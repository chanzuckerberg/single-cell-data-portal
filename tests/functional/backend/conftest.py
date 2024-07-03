import os

import pytest
import requests
from functional.backend.utils import make_dp_auth_header, upload_and_wait
from requests.adapters import HTTPAdapter, Retry

from backend.common.corpora_config import CorporaAuthConfig
from tests.functional.backend.constants import API_URL
from tests.functional.backend.utils import distributed_singleton, get_auth_token, make_cookie


@pytest.fixture(scope="session")
def config():
    return CorporaAuthConfig()


@pytest.fixture(scope="session")
def deployment_stage():
    return os.environ["DEPLOYMENT_STAGE"]


@pytest.fixture(scope="session")
def dataset_uri():
    return (
        "https://www.dropbox.com/scl/fi/y50umqlcrbz21a6jgu99z/5_0_0_example_valid.h5ad?rlkey"
        "=s7p6ybyx082hswix26hbl11pm&dl=0"
    )


@pytest.fixture(scope="session")
def proxy_auth_token(config, deployment_stage, tmp_path_factory, worker_id) -> dict:
    """
    Generate a proxy token for rdev. If running in parallel mode this will be shared across workers to avoid rate
    limiting
    """

    def _proxy_auth_token() -> dict:
        if deployment_stage == "rdev":
            payload = {
                "client_id": config.test_app_id,
                "client_secret": config.test_app_secret,
                "grant_type": "client_credentials",
                "audience": "https://api.cellxgene.dev.single-cell.czi.technology/dp/v1/curator",
            }
            headers = {"content-type": "application/json"}
            res = requests.post("https://czi-cellxgene-dev.us.auth0.com/oauth/token", json=payload, headers=headers)
            res.raise_for_status()
            access_token = res.json()["access_token"]
            return {"Authorization": f"Bearer {access_token}"}
        return {}

    return distributed_singleton(tmp_path_factory, worker_id, _proxy_auth_token)


@pytest.fixture(scope="session")
def session(proxy_auth_token):
    session = requests.Session()
    retry_config = Retry(
        total=7,
        backoff_factor=2,
        status_forcelist=[500, 502, 503, 504],
        allowed_methods={"DELETE", "GET", "HEAD", "PUT", "POST"},
    )
    session.mount("https://", HTTPAdapter(max_retries=retry_config))
    session.headers.update(**proxy_auth_token)
    yield session
    session.close()


@pytest.fixture(scope="session")
def functest_auth_token(config, session, deployment_stage, tmp_path_factory, worker_id):
    def _auth_token():
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
def upload_dataset(session, api_url, curator_cookie, deployment_stage, request):
    def _upload_dataset(collection_id, dropbox_url, existing_dataset_id=None, cleanup=True):
        result = upload_and_wait(
            session, api_url, curator_cookie, deployment_stage, collection_id, dropbox_url, existing_dataset_id
        )
        dataset_id = result["dataset_id"]
        if cleanup:
            request.addfinalizer(
                lambda: session.delete(
                    f"{api_url}/dp/v1/datasets/{dataset_id}", headers=make_dp_auth_header(curator_cookie)
                )
            )
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
