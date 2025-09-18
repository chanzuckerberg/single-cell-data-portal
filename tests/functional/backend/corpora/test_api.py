import json

import pytest
import requests
from requests import HTTPError

from backend.common.constants import DATA_SUBMISSION_POLICY_VERSION
from tests.functional.backend.constants import ATAC_SEQ_MANIFEST, DATASET_MANIFEST, DATASET_URI, VISIUM_DATASET_URI
from tests.functional.backend.skip_reason import skip_creation_on_prod
from tests.functional.backend.utils import assertStatusCode, create_test_collection


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
