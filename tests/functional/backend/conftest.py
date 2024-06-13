import json
import os
import time

import pytest
import requests
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
def visium_dataset_uri():
    return (
        "https://www.dropbox.com/scl/fi/lmhue0va6ihk50ivp26da/visium_small.h5ad?rlkey=n0fo4dyi1ah7ckg9kgzwlhm8s&st"
        "=p7jzej8j&dl=0"
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
    def _curation_api_access_token(session, api_url, config):
        response = session.post(
            f"{api_url}/curation/v1/auth/token",
            headers={"x-api-key": config.super_curator_api_key},
        )
        response.raise_for_status()
        return response.json()["access_token"]

    return distributed_singleton(tmp_path_factory, worker_id, _curation_api_access_token)


@pytest.fixture(scope="session")
def upload_and_wait(session, api_url, curator_cookie, deployment_stage, request):
    def _upload_and_wait(collection_id, dropbox_url, existing_dataset_id=None, cleanup=True):
        headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
        body = {"url": dropbox_url}

        if existing_dataset_id is None:
            res = session.post(
                f"{api_url}/dp/v1/collections/{collection_id}/upload-links", data=json.dumps(body), headers=headers
            )
        else:
            body["id"] = existing_dataset_id
            res = session.put(
                f"{api_url}/dp/v1/collections/{collection_id}/upload-links", data=json.dumps(body), headers=headers
            )

        res.raise_for_status()
        dataset_id = json.loads(res.content)["dataset_id"]
        if cleanup:
            request.addfinalizer(lambda: session.delete(f"{api_url}/dp/v1/datasets/{dataset_id}", headers=headers))
        assert res.status_code == requests.codes.accepted

        keep_trying = True
        expected_upload_statuses = ["WAITING", "UPLOADING", "UPLOADED"]
        expected_conversion_statuses = ["CONVERTING", "CONVERTED", "FAILED", "UPLOADING", "UPLOADED", "NA", None]
        timer = time.time()
        while keep_trying:
            res = session.get(f"{api_url}/dp/v1/datasets/{dataset_id}/status", headers=headers)
            res.raise_for_status()
            data = json.loads(res.content)
            upload_status = data["upload_status"]
            if upload_status:
                assert upload_status in expected_upload_statuses

            if upload_status == "UPLOADED":
                cxg_status = data.get("cxg_status")
                rds_status = data.get("rds_status")
                h5ad_status = data.get("h5ad_status")
                assert data.get("cxg_status") in expected_conversion_statuses
                if cxg_status == "FAILED":
                    pytest.fail(f"CXG CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
                if rds_status == "FAILED":
                    pytest.fail(f"RDS CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
                if h5ad_status == "FAILED":
                    pytest.fail(f"Anndata CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
                if cxg_status == rds_status == h5ad_status == "UPLOADED":
                    keep_trying = False
            if time.time() >= timer + 1200:
                raise TimeoutError(
                    f"Dataset upload or conversion timed out after 10 min. Check logs for dataset: {dataset_id}"
                )
            time.sleep(10)
        return dataset_id

    return _upload_and_wait


def generate_collection_data(request):
    return {
        "contact_email": "lisbon@gmail.com",
        "contact_name": "Madrid Sparkle",
        "curator_name": "John Smith",
        "description": "Well here are some words",
        "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "https://protocol.com"}],
        "name": request.function.__name__,
    }
