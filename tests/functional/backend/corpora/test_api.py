import json

import pytest
import requests
from requests import HTTPError

from backend.common.constants import DATA_SUBMISSION_POLICY_VERSION
from tests.functional.backend.constants import (
    ATAC_SEQ_MANIFEST,
    DATASET_MANIFEST,
    DATASET_URI,
    PRE_ANALYSIS_DATASET_MANIFEST,
    VISIUM_DATASET_URI,
)
from tests.functional.backend.skip_reason import skip_creation_on_prod
from tests.functional.backend.utils import (
    assertStatusCode,
    create_test_collection,
    wait_for_dataset_processing_via_curation_api,
)


@skip_creation_on_prod
def test_version(session, api_url):
    res = session.get(f"{api_url}/dp/v1/deployed_version")
    res.raise_for_status()
    assert res.status_code == requests.codes.ok
    assert len(res.json()["Data Portal"]) > 0


def test_auth(session, api_url, curator_cookie):
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/userinfo", headers=headers)
    res.raise_for_status()
    assert res.status_code == requests.codes.ok
    data = json.loads(res.content)
    assert data["email"] == "functest@example.com"


def test_root_route(session, api_url):
    res = session.get(f"{api_url}/")
    res.raise_for_status()
    assert res.status_code == requests.codes.ok


def test_get_collections(session, api_url):
    res = session.get(f"{api_url}/dp/v1/collections")
    res.raise_for_status()
    assert res.status_code == requests.codes.ok
    data = json.loads(res.content)
    for collection in data["collections"]:
        assert isinstance(collection["id"], str)
        assert isinstance(collection["created_at"], float)


@skip_creation_on_prod
def test_collection_flow(
    session,
    api_url,
    curator_cookie,
    upload_dataset,
    upload_collection_metadata,
    collection_data_DOI_update,
    collection_data,
):
    # create collection
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    res = session.post(f"{api_url}/dp/v1/collections", data=json.dumps(collection_data), headers=headers)
    res.raise_for_status()
    data = json.loads(res.content)
    collection_id = data["collection_id"]
    assertStatusCode(requests.codes.created, res)
    assert "collection_id" in data

    # Test created collection is private
    res = session.get(f"{api_url}/dp/v1/collections", headers=headers)
    data = json.loads(res.content)
    private_collection_ids = []
    for collection in data["collections"]:
        if collection["visibility"] == "PRIVATE":
            private_collection_ids.append(collection["id"])
    assert collection_id in private_collection_ids

    # Test update collection info
    updated_data = {
        "contact_email": "person@random.com",
        "contact_name": "Doctor Who",
        "description": "These are different words",
        "links": [{"link_name": "The Source", "link_type": "DATA_SOURCE", "link_url": "https://datasource.com"}],
        "name": "lots of cells",
    }
    res = session.put(f"{api_url}/dp/v1/collections/{collection_id}", data=json.dumps(updated_data), headers=headers)
    res.raise_for_status()
    data = json.loads(res.content)
    data.pop("access_type")
    for key in updated_data:
        assert updated_data[key] == data[key]

    upload_dataset(collection_id, DATASET_URI)

    # Test collection DOI update and ensure dataset updates are triggered
    upload_collection_metadata(collection_id, collection_data_DOI_update)

    # make collection public
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(f"{api_url}/dp/v1/collections/{collection_id}/publish", headers=headers, data=json.dumps(body))
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # get canonical collection_id
    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    data = json.loads(res.content)
    canonical_collection_id = data["id"]

    # check collection returns as public
    res = session.get(f"{api_url}/dp/v1/collections", headers=headers)
    data = json.loads(res.content)
    public_collection_ids = []
    for collection in data["collections"]:
        if collection["visibility"] == "PUBLIC":
            public_collection_ids.append(collection["id"])

    assert canonical_collection_id in public_collection_ids

    # Test everyone can retrieve a public collection
    no_auth_headers = {"Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/collections", headers=no_auth_headers)
    data = json.loads(res.content)
    collection_ids = [x["id"] for x in data["collections"]]
    assert canonical_collection_id in collection_ids

    # Test a public collection cannot be tombstoned
    res = session.delete(f"{api_url}/dp/v1/collections/{canonical_collection_id}", headers=headers)
    assertStatusCode(requests.codes.method_not_allowed, res)

    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    assertStatusCode(requests.codes.ok, res)


@skip_creation_on_prod
def test_delete_private_collection(session, api_url, curator_cookie, collection_data, request):
    # create collection
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    res = session.post(f"{api_url}/dp/v1/collections", data=json.dumps(collection_data), headers=headers)
    res.raise_for_status()
    data = json.loads(res.content)
    collection_id = data["collection_id"]
    request.addfinalizer(lambda: session.delete(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers))
    assertStatusCode(requests.codes.created, res)
    assert "collection_id" in data

    # check created collection returns as private
    res = session.get(f"{api_url}/dp/v1/collections", headers=headers)
    data = json.loads(res.content)
    private_collection_ids = []
    for collection in data["collections"]:
        if collection["visibility"] == "PRIVATE":
            private_collection_ids.append(collection["id"])

    assert collection_id in private_collection_ids

    # delete collection
    res = session.delete(f"{api_url}/dp/v1/collections/{collection_id}?visibility=PRIVATE", headers=headers)
    res.raise_for_status()
    assertStatusCode(requests.codes.no_content, res)

    # check collection gone
    no_auth_headers = {"Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/collections?visibility=PRIVATE", headers=no_auth_headers)
    data = json.loads(res.content)
    collection_ids = [x["id"] for x in data["collections"]]
    assert collection_id not in collection_ids


@skip_creation_on_prod
def test_dataset_upload_flow_with_dataset(
    session,
    curator_cookie,
    api_url,
    upload_dataset,
    upload_dataset_title,
    request,
    collection_data,
    dataset_title_update,
):
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    collection_id = create_test_collection(headers, request, session, api_url, collection_data)
    _verify_upload_and_delete_succeeded(
        collection_id,
        headers,
        DATASET_URI,
        dataset_title_update,
        session,
        api_url,
        upload_dataset,
        upload_dataset_title,
    )


@skip_creation_on_prod
def test_dataset_upload_flow_with_visium_dataset(
    session,
    curator_cookie,
    api_url,
    upload_dataset,
    upload_dataset_title,
    request,
    collection_data,
    dataset_title_update,
):
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    collection_id = create_test_collection(headers, request, session, api_url, collection_data)
    _verify_upload_and_delete_succeeded(
        collection_id,
        headers,
        VISIUM_DATASET_URI,
        dataset_title_update,
        session,
        api_url,
        upload_dataset,
        upload_dataset_title,
    )


@skip_creation_on_prod
def test_dataset_upload_flow_with_atac_seq_dataset(
    session,
    curator_cookie,
    api_url,
    upload_manifest,
    upload_dataset_title,
    request,
    collection_data,
    dataset_title_update,
    curation_api_access_token,
):
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    collection_id = create_test_collection(
        headers,
        request,
        session,
        api_url,
        collection_data,
    )
    _verify_upload_and_delete_succeeded(
        collection_id,
        headers,
        ATAC_SEQ_MANIFEST,
        dataset_title_update,
        session,
        api_url,
        upload_manifest,
        upload_dataset_title,
    )


@skip_creation_on_prod
def test_dataset_reupload_flow_from_manifest(
    session,
    curator_cookie,
    api_url,
    upload_manifest,
    upload_dataset_title,
    request,
    collection_data,
    dataset_title_update,
    curation_api_access_token,
):
    """Test reupload from public urls."""
    headers_dp = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    headers_curation = {"Authorization": f"Bearer {curation_api_access_token}", "Content-Type": "application/json"}
    collection_id = create_test_collection(
        headers_dp,
        request,
        session,
        api_url,
        collection_data,
    )
    result = upload_manifest(collection_id, DATASET_MANIFEST)
    dataset_id = result["dataset_id"]

    # publish the collection
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(
        f"{api_url}/dp/v1/collections/{collection_id}/publish", headers=headers_dp, data=json.dumps(body)
    )
    assertStatusCode(requests.codes.accepted, res)

    # start a revision
    res = session.post(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers_dp)
    assertStatusCode(201, res)
    collection_id = res.json()["id"]

    # get the manifest and ensure it has expected content
    resp = session.get(
        f"{api_url}/curation/v1/collections/{collection_id}/datasets/{dataset_id}/manifest",
        headers=headers_curation,
    )
    assertStatusCode(200, resp)
    new_manifest = resp.json()
    assert set(new_manifest.keys()) == set(
        DATASET_MANIFEST.keys()
    ), f"Manifest keys do not match expected, {new_manifest=}"
    # re-upload the manifest from the public urls to ensure re-upload works as expected
    upload_manifest(collection_id, new_manifest, existing_dataset_id=dataset_id)

    # publish the revision
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(
        f"{api_url}/dp/v1/collections/{collection_id}/publish", headers=headers_dp, data=json.dumps(body)
    )
    assertStatusCode(requests.codes.accepted, res)


def _verify_upload_and_delete_succeeded(
    collection_id,
    headers,
    req_body,
    title_update_body,
    session,
    api_url,
    upload_and_wait,
    upload_dataset_title_and_wait,
):
    result = upload_and_wait(collection_id, req_body)
    dataset_id = result["dataset_id"]
    version_id = result["version_id"]
    # test non owner cant retrieve status
    no_auth_headers = {"Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/datasets/{dataset_id}/status", headers=no_auth_headers)
    with pytest.raises(HTTPError):
        res.raise_for_status()

    # update title and await dataset update
    updated_dataset_id = upload_dataset_title_and_wait(collection_id, version_id, title_update_body)

    # Test dataset deletion
    res = session.delete(f"{api_url}/dp/v1/datasets/{updated_dataset_id}", headers=headers)
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # Check that the dataset is gone from collection version
    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    data = json.loads(res.content)
    datasets = data["datasets"]
    dataset_ids = [dataset.get("id") for dataset in datasets]
    assert updated_dataset_id not in dataset_ids


# --- is_pre_analysis functional tests ---


@skip_creation_on_prod
def test_pre_analysis_collection_create_and_get(session, api_url, curation_api_access_token, request):
    """Creating a collection with is_pre_analysis=True via the Curation API stores and returns the flag."""
    headers = {"Authorization": f"Bearer {curation_api_access_token}", "Content-Type": "application/json"}
    body = {
        "contact_email": "functest@example.com",
        "contact_name": "Func Test",
        "description": "pre-analysis functional test",
        "name": "test_pre_analysis_collection_create_and_get",
        "is_pre_analysis": True,
    }
    res = session.post(f"{api_url}/curation/v1/collections", data=json.dumps(body), headers=headers)
    assertStatusCode(requests.codes.created, res)
    collection_id = res.json()["collection_id"]
    request.addfinalizer(lambda: session.delete(f"{api_url}/curation/v1/collections/{collection_id}", headers=headers))

    res = session.get(f"{api_url}/curation/v1/collections/{collection_id}", headers=headers)
    assertStatusCode(requests.codes.ok, res)
    data = res.json()
    assert "is_pre_analysis" in data
    assert data["is_pre_analysis"] is True


@skip_creation_on_prod
def test_pre_analysis_defaults_false(session, api_url, curation_api_access_token, request):
    """Creating a collection without is_pre_analysis returns is_pre_analysis=False."""
    headers = {"Authorization": f"Bearer {curation_api_access_token}", "Content-Type": "application/json"}
    body = {
        "contact_email": "functest@example.com",
        "contact_name": "Func Test",
        "description": "pre-analysis default test",
        "name": "test_pre_analysis_defaults_false",
    }
    res = session.post(f"{api_url}/curation/v1/collections", data=json.dumps(body), headers=headers)
    assertStatusCode(requests.codes.created, res)
    collection_id = res.json()["collection_id"]
    request.addfinalizer(lambda: session.delete(f"{api_url}/curation/v1/collections/{collection_id}", headers=headers))

    res = session.get(f"{api_url}/curation/v1/collections/{collection_id}", headers=headers)
    assertStatusCode(requests.codes.ok, res)
    data = res.json()
    assert "is_pre_analysis" in data
    assert data["is_pre_analysis"] is False


@skip_creation_on_prod
def test_pre_analysis_patch_rejected(session, api_url, curation_api_access_token, request):
    """Attempting to PATCH is_pre_analysis on an existing collection returns 405."""
    headers = {"Authorization": f"Bearer {curation_api_access_token}", "Content-Type": "application/json"}
    body = {
        "contact_email": "functest@example.com",
        "contact_name": "Func Test",
        "description": "patch rejected test",
        "name": "test_pre_analysis_patch_rejected",
    }
    res = session.post(f"{api_url}/curation/v1/collections", data=json.dumps(body), headers=headers)
    assertStatusCode(requests.codes.created, res)
    collection_id = res.json()["collection_id"]
    request.addfinalizer(lambda: session.delete(f"{api_url}/curation/v1/collections/{collection_id}", headers=headers))

    res = session.patch(
        f"{api_url}/curation/v1/collections/{collection_id}",
        data=json.dumps({"is_pre_analysis": True}),
        headers=headers,
    )
    assertStatusCode(405, res)


def test_get_datasets_analysis_filter(session, api_url):
    """GET /curation/v1/datasets with ?analysis= returns only matching datasets."""
    # pre-analysis filter — all returned datasets must have is_pre_analysis=True
    res = session.get(f"{api_url}/curation/v1/datasets?analysis=pre-analysis")
    assertStatusCode(requests.codes.ok, res)
    for dataset in res.json():
        assert dataset.get("is_pre_analysis") is True

    # post-analysis filter — all returned datasets must have is_pre_analysis=False
    res = session.get(f"{api_url}/curation/v1/datasets?analysis=post-analysis")
    assertStatusCode(requests.codes.ok, res)
    for dataset in res.json():
        assert dataset.get("is_pre_analysis") is False

    # invalid value — expect 400
    res = session.get(f"{api_url}/curation/v1/datasets?analysis=invalid")
    assertStatusCode(400, res)


def test_get_datasets_includes_is_pre_analysis_field(session, api_url):
    """Every dataset returned by GET /curation/v1/datasets includes is_pre_analysis."""
    res = session.get(f"{api_url}/curation/v1/datasets")
    assertStatusCode(requests.codes.ok, res)
    for dataset in res.json():
        assert "is_pre_analysis" in dataset


# --- pre-analysis upload functional tests ---


@skip_creation_on_prod
def test_pre_analysis_dataset_upload_and_process(
    session,
    api_url,
    curation_api_access_token,
    curator_cookie,
    request,
):
    """
    Upload a pre-analysis h5ad to a pre-analysis collection and verify the full
    pipeline completes (CXG conversion is skipped) and all API responses reflect
    is_pre_analysis=True.

    Everything uses the curation Bearer token — the same path a real API user takes.
    Status is polled via GET /curation/v1/collections/{id} (processing_status field)
    rather than the dp/v1 status endpoint, which is cookie-only and unavailable to
    Bearer-token callers.
    """
    headers = {"Authorization": f"Bearer {curation_api_access_token}", "Content-Type": "application/json"}

    # Create a pre-analysis collection.
    res = session.post(
        f"{api_url}/curation/v1/collections",
        data=json.dumps(
            {
                "contact_email": "functest@example.com",
                "contact_name": "Func Test",
                "description": "pre-analysis upload functional test",
                "name": "test_pre_analysis_dataset_upload_and_process",
                "is_pre_analysis": True,
            }
        ),
        headers=headers,
    )
    assertStatusCode(requests.codes.created, res)
    collection_id = res.json()["collection_id"]
    request.addfinalizer(lambda: session.delete(f"{api_url}/curation/v1/collections/{collection_id}", headers=headers))

    # Create a dataset slot and submit the manifest.
    res = session.post(f"{api_url}/curation/v1/collections/{collection_id}/datasets", headers=headers)
    assertStatusCode(201, res)
    dataset_id = res.json()["dataset_id"]

    res = session.put(
        f"{api_url}/curation/v1/collections/{collection_id}/datasets/{dataset_id}/manifest",
        data=json.dumps(PRE_ANALYSIS_DATASET_MANIFEST),
        headers=headers,
    )
    assertStatusCode(202, res)

    # Wait for processing via the curation API (Bearer-token-accessible).
    result = wait_for_dataset_processing_via_curation_api(
        session, api_url, curation_api_access_token, collection_id, dataset_id
    )
    assert not result["errors"], f"Pre-analysis dataset processing failed: {result['errors']}"

    # Verify is_pre_analysis on the collection and dataset responses.
    res = session.get(f"{api_url}/curation/v1/collections/{collection_id}", headers=headers)
    assertStatusCode(200, res)
    collection_data = res.json()
    assert collection_data.get("is_pre_analysis") is True
    dataset_entry = next(d for d in collection_data["datasets"] if d["dataset_id"] == dataset_id)
    assert dataset_entry.get("is_pre_analysis") is True

    res = session.get(
        f"{api_url}/curation/v1/collections/{collection_id}/datasets/{dataset_id}",
        headers=headers,
    )
    assertStatusCode(200, res)
    assert res.json().get("is_pre_analysis") is True

    # Verify the dataset index filters work correctly.
    # The collection is still private (unpublished), so visibility=PRIVATE is required.
    res = session.get(f"{api_url}/curation/v1/datasets?analysis=pre-analysis&visibility=PRIVATE", headers=headers)
    assertStatusCode(200, res)
    assert dataset_id in [d["dataset_id"] for d in res.json()]

    res = session.get(f"{api_url}/curation/v1/datasets?analysis=post-analysis&visibility=PRIVATE", headers=headers)
    assertStatusCode(200, res)
    assert dataset_id not in [d["dataset_id"] for d in res.json()]

    # Publish the collection, then verify the public analysis filter works without visibility=PRIVATE.
    headers_dp = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(
        f"{api_url}/dp/v1/collections/{collection_id}/publish", headers=headers_dp, data=json.dumps(body)
    )
    assertStatusCode(requests.codes.accepted, res)

    # After publishing the collection is public, so the filter should work without visibility=PRIVATE.
    res = session.get(f"{api_url}/curation/v1/datasets?analysis=pre-analysis", headers=headers)
    assertStatusCode(200, res)
    assert dataset_id in [d["dataset_id"] for d in res.json()]

    res = session.get(f"{api_url}/curation/v1/datasets?analysis=post-analysis", headers=headers)
    assertStatusCode(200, res)
    assert dataset_id not in [d["dataset_id"] for d in res.json()]
