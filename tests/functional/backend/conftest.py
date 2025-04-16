import os

import pytest

from backend.common.corpora_config import CorporaAuthConfig
from tests.functional.backend.constants import API_URL
from tests.functional.backend.distributed import distributed_singleton
from tests.functional.backend.utils import (
    get_auth_token,
    get_curation_api_access_token,
    make_cookie,
    make_proxy_auth_token,
    make_session,
    update_metadata_and_wait,
    update_title_and_wait,
    upload_manifest_and_wait,
    upload_url_and_wait,
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
    def _functest_auth_token() -> dict[str, str]:
        username = config.functest_account_username
        password = config.functest_account_password
        return get_auth_token(username, password, session, config, deployment_stage)

    return distributed_singleton(tmp_path_factory, worker_id, _functest_auth_token)


@pytest.fixture(scope="session")
def curator_cookie(functest_auth_token):
    return make_cookie(functest_auth_token)


@pytest.fixture(scope="session")
def api_url(deployment_stage):
    return API_URL.get(deployment_stage)


@pytest.fixture(scope="session")
def curation_api_access_token(session, api_url, config, tmp_path_factory, worker_id):
    def _curation_api_access_token() -> str:
        return get_curation_api_access_token(session, api_url, config)

    return distributed_singleton(tmp_path_factory, worker_id, _curation_api_access_token)


@pytest.fixture(scope="session")
def upload_dataset(session, api_url, curator_cookie, request):
    def _upload_dataset(collection_id, dropbox_url, existing_dataset_id=None):
        result = upload_url_and_wait(session, api_url, curator_cookie, collection_id, dropbox_url, existing_dataset_id)
        dataset_id = result["dataset_id"]
        headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
        request.addfinalizer(lambda: session.delete(f"{api_url}/dp/v1/datasets/{dataset_id}", headers=headers))
        if result["errors"]:
            raise pytest.fail(str(result["errors"]))
        return dataset_id

    return _upload_dataset


@pytest.fixture(scope="session")
def upload_collection_metadata(session, api_url, curator_cookie, request):
    def _upload_collection_metadata(collection_id, metadata):
        collection_errors = update_metadata_and_wait(session, api_url, curator_cookie, collection_id, metadata)
        if any(errors for errors in collection_errors.values()):
            raise pytest.fail(str(collection_errors))
        dataset_ids = list(collection_errors.keys())
        return collection_id, dataset_ids

    return _upload_collection_metadata


@pytest.fixture(scope="session")
def upload_dataset_title(session, api_url, curator_cookie, request):
    def _upload_dataset_title(collection_id, title):
        return update_title_and_wait(session, api_url, curator_cookie, collection_id, title)

    return _upload_dataset_title


@pytest.fixture(scope="session")
def upload_manifest(session, api_url, curation_api_access_token, curator_cookie, request):
    def _upload_manifest(collection_id: str, manifest: dict, existing_dataset_id=None):
        result = upload_manifest_and_wait(
            session, api_url, curation_api_access_token, curator_cookie, collection_id, manifest, existing_dataset_id
        )
        dataset_id = result["dataset_id"]
        headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
        request.addfinalizer(lambda: session.delete(f"{api_url}/dp/v1/datasets/{dataset_id}", headers=headers))
        if result["errors"]:
            raise pytest.fail(str(result["errors"]))
        return dataset_id

    return _upload_manifest


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


@pytest.fixture()
def collection_data_DOI_update(request):
    return {
        "contact_email": "lisbon@gmail.com",
        "contact_name": "Madrid Sparkle",
        "curator_name": "John Smith",
        "description": "Well here are some words",
        "links": [
            {"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "https://protocol.com"},
            {"link_name": "", "link_type": "DOI", "link_url": "10.1093/nar/gkae1142"},
        ],
        "name": request.function.__name__,
    }


@pytest.fixture()
def dataset_title_update(request):
    return {"title": f"Updated Title for {request.function.__name__}"}
