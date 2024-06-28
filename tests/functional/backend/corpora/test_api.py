import json
import os

import pytest
import requests
from requests import HTTPError

from backend.common.constants import DATA_SUBMISSION_POLICY_VERSION
from tests.functional.backend.utils import assertStatusCode


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


@pytest.mark.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
def test_collection_flow(session, api_url, curator_cookie, upload_and_wait, dataset_uri, collection_data):
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

    upload_and_wait(collection_id, dataset_uri)

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


@pytest.mark.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
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


@pytest.mark.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
def test_dataset_upload_flow_with_visium_dataset(
    session, curator_cookie, api_url, upload_and_wait, visium_dataset_uri, request
):
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    collection_id = _create_test_collection(
        headers, request, session, api_url, "test_dataset_upload_flow_with_visium_dataset"
    )
    _verify_upload_and_delete_succeeded(collection_id, headers, visium_dataset_uri, session, api_url, upload_and_wait)


def _create_test_collection(headers, request, session, api_url, name="my2collection"):
    body = {
        "contact_email": "lisbon@gmail.com",
        "contact_name": "Madrid Sparkle",
        "curator_name": "John Smith",
        "description": "Well here are some words",
        "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "https://protocol.com"}],
        "name": name,
    }

    res = session.post(f"{api_url}/dp/v1/collections", data=json.dumps(body), headers=headers)
    res.raise_for_status()
    data = json.loads(res.content)
    collection_id = data["collection_id"]
    request.addfinalizer(lambda: session.delete(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers))
    assertStatusCode(requests.codes.created, res)
    return collection_id


def _verify_upload_and_delete_succeeded(collection_id, headers, dataset_uri, session, api_url, upload_and_wait):
    dataset_id = upload_and_wait(collection_id, dataset_uri)
    # test non owner cant retrieve status
    no_auth_headers = {"Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/datasets/{dataset_id}/status", headers=no_auth_headers)
    with pytest.raises(HTTPError):
        res.raise_for_status()

    # Test dataset deletion
    res = session.delete(f"{api_url}/dp/v1/datasets/{dataset_id}", headers=headers)
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # Check that the dataset is gone from collection version
    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    data = json.loads(res.content)
    datasets = data["datasets"]
    dataset_ids = [dataset.get("id") for dataset in datasets]
    assert dataset_id not in dataset_ids
