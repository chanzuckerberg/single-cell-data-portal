import base64
import json
import time
from typing import Optional

import requests
from requests import Session
from requests.adapters import HTTPAdapter, Retry

from backend.common.corpora_config import CorporaAuthConfig
from tests.functional.backend.constants import AUDIENCE


def get_auth_token(
    username: str,
    password: str,
    session: Session,
    config: CorporaAuthConfig,
    deployment_stage: str,
    additional_claims: Optional[list] = None,
) -> dict[str, str]:
    standard_claims = "openid profile email offline"
    if additional_claims:
        additional_claims.append(standard_claims)
        claims = " ".join(additional_claims)
    else:
        claims = standard_claims
    response = session.post(
        "https://czi-cellxgene-dev.us.auth0.com/oauth/token",  # hardcoded becasue this is only needed for dev and rdev
        headers={"content-type": "application/x-www-form-urlencoded"},
        data=dict(
            grant_type="password",
            username=username,
            password=password,
            audience=AUDIENCE.get(deployment_stage),
            scope=claims,
            client_id=config.client_id,
            client_secret=config.client_secret,
        ),
    )
    response.raise_for_status()
    access_token = response.json()["access_token"]
    id_token = response.json()["id_token"]
    token = {"access_token": access_token, "id_token": id_token}
    return token


def make_cookie(auth_token: dict) -> str:
    return base64.b64encode(json.dumps(auth_token).encode("utf-8")).decode()


def assertStatusCode(expected_status: int, actual_resp: requests.Response):
    request_id = actual_resp.headers.get("X-Request-Id")
    assert expected_status == actual_resp.status_code, f"{request_id=}, {actual_resp.content=}"


def create_test_collection(headers, request, session, api_url, body):
    res = session.post(f"{api_url}/dp/v1/collections", data=json.dumps(body), headers=headers)
    assertStatusCode(requests.codes.created, res)
    data = json.loads(res.content)
    collection_id = data["collection_id"]
    request.addfinalizer(lambda: session.delete(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers))
    return collection_id


def create_explorer_url(dataset_id: str, deployment_stage: str) -> str:
    return f"https://cellxgene.{deployment_stage}.single-cell.czi.technology/e/{dataset_id}.cxg/"


def update_metadata_and_wait(session, api_url, curator_cookie, collection_id, metadata):
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    res = session.put(f"{api_url}/dp/v1/collections/{collection_id}", data=json.dumps(metadata), headers=headers)
    res.raise_for_status()
    # ensure metadata update is queued for each dataset
    collection = json.loads(res.content)
    dataset_ids = [dataset["id"] for dataset in collection["datasets"]]
    for dataset_id in dataset_ids:
        res = session.get(f"{api_url}/dp/v1/datasets/{dataset_id}/status", headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        assert data["processing_status"] == "INITIALIZED"

    collection_errors = {}
    for dataset_id in dataset_ids:
        result = _wait_for_dataset_status(session, api_url, dataset_id, headers)
        dataset_id = result["dataset_id"]
        dataset_errors = result["errors"]
        collection_errors[dataset_id] = dataset_errors

    return collection_errors


def update_title_and_wait(session, api_url, curator_cookie, collection_id, dataset_id, dataset_title_update):
    headers = {
        "Cookie": f"cxguser={curator_cookie}",
        "Content-Type": "application/json",
    }

    patch_url = f"{api_url}/dp/v1/collections/{collection_id}/datasets/{dataset_id}"
    res = session.patch(patch_url, data=json.dumps(dataset_title_update), headers=headers)
    res.raise_for_status()

    # ensure metadata update is queued for dataset
    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    res.raise_for_status()
    collection = json.loads(res.content)
    updated_dataset_id = [dataset["id"] for dataset in collection["datasets"]][0]
    res = session.get(f"{api_url}/dp/v1/datasets/{updated_dataset_id}/status", headers=headers)
    res.raise_for_status()
    data = json.loads(res.content)
    assert data["processing_status"] == "INITIALIZED"

    result = _wait_for_dataset_status(session, api_url, updated_dataset_id, headers)
    dataset_errors = result["errors"]

    # Check if title was updated
    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    res.raise_for_status()
    collection = json.loads(res.content)
    updated_dataset = next((dataset for dataset in collection["datasets"] if dataset["id"] == updated_dataset_id), None)
    assert updated_dataset["name"] == dataset_title_update["title"]
    assert not dataset_errors

    return updated_dataset_id


def upload_url_and_wait(session, api_url, curator_cookie, collection_id, dropbox_url, existing_dataset_id=None):
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
    assertStatusCode(requests.codes.accepted, res)
    dataset_id = json.loads(res.content)["dataset_id"]

    return _wait_for_dataset_status(session, api_url, dataset_id, headers)


def upload_manifest_and_wait(
    session, api_url, curation_api_access_token, curator_cookie, collection_id, manifest, existing_dataset_id=None
):
    headers = {"Authorization": f"Bearer {curation_api_access_token}", "Content-Type": "application/json"}

    if not existing_dataset_id:
        # Create dataset id
        res = session.post(f"{api_url}/curation/v1/collections/{collection_id}/datasets", headers=headers)
        assertStatusCode(201, res)
        dataset_id = json.loads(res.content)["dataset_id"]
        res = session.get(f"{api_url}/curation/v1/collections/{collection_id}", headers=headers)
        assertStatusCode(200, res)
        version_id = json.loads(res.content)["datasets"][0]["dataset_version_id"]
    else:
        dataset_id = existing_dataset_id

    # Upload manifest
    res = session.put(
        f"{api_url}/curation/v1/collections/{collection_id}/datasets/{dataset_id}/manifest",
        data=json.dumps(manifest),
        headers=headers,
    )
    assertStatusCode(202, res)

    # Wait for dataset status
    return _wait_for_dataset_status(
        session, api_url, version_id, {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    )


def _wait_for_dataset_status(session, api_url, dataset_id, headers):
    errors = []
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
            atac_status = data.get("atac_status")
            assert data.get("cxg_status") in expected_conversion_statuses
            if cxg_status == "FAILED":
                errors.append(f"CXG CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
            if h5ad_status == "FAILED":
                errors.append(f"Anndata CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
            if atac_status == "FAILED":
                errors.append(f"Atac CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
            if rds_status == "FAILED":
                errors.append(f"RDS CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
            if any(
                [
                    (cxg_status == h5ad_status == "UPLOADED" or cxg_status == h5ad_status == "CONVERTED")
                    and rds_status == "SKIPPED"
                    and atac_status in ["SKIPPED", "UPLOADED", "NA", "COPIED"],
                    errors,
                ]
            ):
                keep_trying = False
            if time.time() >= timer + 3600:
                raise TimeoutError(
                    f"Dataset upload or conversion timed out after 30 min. Check logs for dataset: {dataset_id}"
                    f"upload_status: {upload_status} cxg_status: {cxg_status}, rds_status: {rds_status}, "
                    f"h5ad_status: {h5ad_status}, atac_status: {atac_status}"
                )
        time.sleep(10)
    return {"dataset_id": dataset_id, "errors": errors}


def make_session(proxy_auth_token):
    retry_strategy = Retry(
        total=7,
        backoff_factor=2,
        status_forcelist=[500, 502, 503, 504],
        allowed_methods={"DELETE", "GET", "HEAD", "PUT", "POST"},
    )
    http_adapter = HTTPAdapter(max_retries=retry_strategy)
    session = requests.Session()
    session.mount("https://", http_adapter)
    session.headers.update(**proxy_auth_token)
    return session


def make_proxy_auth_token(config, deployment_stage) -> dict:
    """
    Generate a proxy token for rdev. If running in parallel mode this will be shared across workers to avoid rate
    limiting
    """

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


def get_curation_api_access_token(session, api_url, config) -> str:
    response = session.post(
        f"{api_url}/curation/v1/auth/token",
        headers={"x-api-key": config.super_curator_api_key},
    )
    response.raise_for_status()
    return response.json()["access_token"]
