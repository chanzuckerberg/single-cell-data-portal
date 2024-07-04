import os

import pytest

from backend.common.corpora_config import CorporaAuthConfig
from tests.functional.backend.constants import API_URL
from tests.functional.backend.distributed import distributed_singleton
from tests.functional.backend.utils import (
    get_auth_token,
    make_cookie,
    make_proxy_auth_token,
    make_session,
    upload_and_wait,
)


@pytest.fixture(scope="session")
def config():
    return CorporaAuthConfig()


@pytest.fixture(scope="session")
def deployment_stage():
    return os.environ["DEPLOYMENT_STAGE"]


@pytest.fixture(scope="session")
def proxy_auth_token(config, deployment_stage, tmp_path_factory, worker_id) -> dict:
    """
    Generate a proxy token for rdev. If running in parallel mode this will be shared across workers to avoid rate
    limiting
    """

    def _proxy_auth_token() -> dict:
        return make_proxy_auth_token(config, deployment_stage)

    return distributed_singleton(tmp_path_factory, worker_id, _proxy_auth_token)


@pytest.fixture(scope="session")
def session(proxy_auth_token):
    session = make_session(proxy_auth_token)
    yield session
    session.close()


@pytest.fixture(scope="session")
def functest_auth_token(config, session, deployment_stage, tmp_path_factory, worker_id):
    def _auth_token() -> dict[str, str]:
        username = config.functest_account_username
        password = config.functest_account_password
        return get_auth_token(username, password, session, config, deployment_stage)

    return distributed_singleton(tmp_path_factory, worker_id, _auth_token)


@pytest.fixture(scope="session")
def curator_cookie(functest_auth_token):
    return make_cookie(functest_auth_token)


@pytest.fixture(scope="session")
def api_url(deployment_stage):
    return API_URL.get(deployment_stage)


@pytest.fixture(scope="session")
def curation_api_access_token(session, api_url, config, tmp_path_factory, worker_id):
    def _curation_api_access_token() -> str:
        response = session.post(
            f"{api_url}/curation/v1/auth/token",
            headers={"x-api-key": config.super_curator_api_key},
        )
        response.raise_for_status()
        return response.json()["access_token"]

    return distributed_singleton(tmp_path_factory, worker_id, _curation_api_access_token)


@pytest.fixture(scope="session")
def upload_dataset(session, api_url, curator_cookie, request):
    def _upload_dataset(collection_id, dropbox_url, existing_dataset_id=None):
        result = upload_and_wait(session, api_url, curator_cookie, collection_id, dropbox_url, existing_dataset_id)
        dataset_id = result["dataset_id"]
        headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
        request.addfinalizer(lambda: session.delete(f"{api_url}/dp/v1/datasets/{dataset_id}", headers=headers))
        if result["errors"]:
            raise pytest.fail(str(result["errors"]))
        return dataset_id

    return _upload_dataset


@pytest.fixture()
def collection_data(request):
    return {
        "contact_email": "lisbon@gmail.com",
        "contact_name": "Madrid Sparkle",
        "curator_name": "John Smith",
        "description": "Well here are some words",
        "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "https://protocol.com"}],
        "name": request.function.__name__,
    }
