import json
import os
import unittest
from urllib.parse import quote

import requests
from tenacity import retry, stop_after_attempt, wait_fixed

from backend.common.constants import DATA_SUBMISSION_POLICY_VERSION
from tests.functional.backend.constants import DATASET_URI
from tests.functional.backend.skip_reason import skip_creation_on_prod, skip_no_explorer_in_rdev
from tests.functional.backend.utils import assertStatusCode, create_explorer_url, create_test_collection


@skip_creation_on_prod
@skip_no_explorer_in_rdev
def test_revision_flow(
    curator_cookie,
    session,
    api_url,
    upload_dataset,
    curation_api_access_token,
    deployment_stage,
    request,
    collection_data,
):
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}

    # create a test collection
    collection_id = create_test_collection(headers, request, session, api_url, collection_data)

    dataset_1_dropbox_url = dataset_2_dropbox_url = DATASET_URI

    # Uploads a dataset
    upload_dataset(collection_id, dataset_1_dropbox_url)

    # make collection public
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(f"{api_url}/dp/v1/collections/{collection_id}/publish", headers=headers, data=json.dumps(body))
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # get canonical collection id, post-publish
    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    data = json.loads(res.content)
    canonical_collection_id = data["id"]

    dataset_response = session.get(f"{api_url}/dp/v1/collections/{canonical_collection_id}").json()["datasets"][0]
    dataset_id = dataset_response["id"]
    explorer_url = dataset_response["dataset_deployments"][0]["url"]

    meta_payload_before_revision_res = session.get(f"{api_url}/dp/v1/datasets/meta?url={explorer_url}")
    meta_payload_before_revision_res.raise_for_status()
    meta_payload_before_revision = meta_payload_before_revision_res.json()

    # Endpoint is eventually consistent
    schema_before_revision = get_schema_with_retries(dataset_id, api_url, session).json()

    # Start a revision
    res = session.post(f"{api_url}/dp/v1/collections/{canonical_collection_id}", headers=headers)
    assertStatusCode(201, res)
    data = json.loads(res.content)
    revision_id = data["id"]

    # "Test updating a dataset in a revision does not effect the published dataset"
    private_dataset_id = res.json()["datasets"][0]["id"]

    meta_payload_res = session.get(f"{api_url}/dp/v1/datasets/meta?url={explorer_url}")
    meta_payload_res.raise_for_status()
    meta_payload = meta_payload_res.json()

    assert meta_payload_before_revision == meta_payload

    # Upload a new dataset
    upload_dataset(
        revision_id,
        dataset_2_dropbox_url,
        existing_dataset_id=private_dataset_id,
    )

    # Check that the published dataset is still the same
    meta_payload_after_revision = session.get(f"{api_url}/dp/v1/datasets/meta?url={explorer_url}").json()
    assert meta_payload_before_revision == meta_payload_after_revision
    schema_after_revision = get_schema_with_retries(dataset_id, api_url, session).json()
    assert schema_before_revision == schema_after_revision

    # Publishing a revised dataset replaces the original dataset
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(f"{api_url}/dp/v1/collections/{revision_id}/publish", headers=headers, data=json.dumps(body))
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    dataset_meta_payload = session.get(f"{api_url}/dp/v1/datasets/meta?url={explorer_url}").json()
    assert dataset_meta_payload["s3_uri"].startswith(f"s3://hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}/")
    assert dataset_meta_payload["s3_uri"].endswith(".cxg/")
    assert (
        dataset_meta_payload["dataset_id"] in dataset_meta_payload["s3_uri"]
    ), "The id of the S3_URI should be the revised dataset id."

    # TODO: add `And the explorer url redirects appropriately`

    # Start a new revision
    res = session.post(f"{api_url}/dp/v1/collections/{canonical_collection_id}", headers=headers)
    assertStatusCode(201, res)
    revision_id = res.json()["id"]

    # Get datasets for the collection (before uploading)
    public_datasets_before = session.get(f"{api_url}/dp/v1/collections/{canonical_collection_id}").json()["datasets"]

    # Upload a new dataset
    another_dataset_id = upload_dataset(revision_id, dataset_1_dropbox_url)

    # Adding a dataset to a revision does not impact public datasets in that collection
    # Get datasets for the collection (after uploading)
    public_datasets_after = session.get(f"{api_url}/dp/v1/collections/{canonical_collection_id}").json()["datasets"]
    unittest.TestCase().assertCountEqual(public_datasets_before, public_datasets_after)

    # Publish the revision
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(f"{api_url}/dp/v1/collections/{revision_id}/publish", headers=headers, data=json.dumps(body))
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # Publishing a revision that contains a new dataset updates the collection page for the data portal (with the new
    # dataset)
    # Check if the last updated dataset_id is among the public datasets
    public_datasets = session.get(f"{api_url}/dp/v1/collections/{canonical_collection_id}").json()["datasets"]
    assert len(public_datasets) == 2
    ids = [dataset["id"] for dataset in public_datasets]
    assert another_dataset_id in ids

    # Start a revision
    res = session.post(f"{api_url}/dp/v1/collections/{canonical_collection_id}", headers=headers)
    assertStatusCode(201, res)
    revision_id = res.json()["id"]

    # This only works if you pick the non replaced dataset.
    dataset_to_delete = res.json()["datasets"][1]
    revision_deleted_dataset_id = dataset_to_delete["id"]
    published_explorer_url = create_explorer_url(revision_deleted_dataset_id, deployment_stage)

    # Delete (tombstone) a dataset (using admin privileges) within the revision
    revision_datasets = session.get(f"{api_url}/curation/v1/collections/{revision_id}").json()["datasets"]
    dataset_id_to_delete = None
    for dataset in revision_datasets:
        if dataset["dataset_version_id"] == revision_deleted_dataset_id:
            dataset_id_to_delete = dataset["dataset_id"]

    curation_api_headers = {"Authorization": f"Bearer {curation_api_access_token}"}
    res = session.delete(
        f"{api_url}/curation/v1/collections/{revision_id}/datasets/{dataset_id_to_delete}?delete_published=true",
        headers=curation_api_headers,
    )
    assertStatusCode(202, res)

    # Deleting a dataset does not effect the published dataset"
    # Check if the dataset is still available
    res = session.get(f"{api_url}/dp/v1/datasets/meta?url={published_explorer_url}")
    assertStatusCode(200, res)

    # Endpoint is eventually consistent
    res = get_schema_with_retries(revision_deleted_dataset_id, api_url, session)
    assertStatusCode(200, res)

    # Publishing a revision that deletes a dataset removes it from the data portal
    # Publish the revision
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(f"{api_url}/dp/v1/collections/{revision_id}/publish", headers=headers, data=json.dumps(body))
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # Check that the dataset doesn't exist anymore
    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    res.raise_for_status()
    datasets = [dataset["id"] for dataset in res.json()["datasets"]]
    assert len(datasets) == 1
    assert revision_deleted_dataset_id not in datasets


def get_schema_with_retries(dataset_id, api_url, session, desired_http_status_code=requests.codes.ok):
    @retry(wait=wait_fixed(2), stop=stop_after_attempt(50))
    def get_s3_uri():
        s3_uri_res = session.get(f"{api_url}/cellxgene/e/{dataset_id}.cxg/api/v0.3/s3_uri", allow_redirects=False)
        assert s3_uri_res.status_code == desired_http_status_code
        return s3_uri_res

    @retry(wait=wait_fixed(2), stop=stop_after_attempt(50))
    def get_schema(s3_uri_response_object):
        # parse s3_uri_response_object content
        s3_path = s3_uri_response_object.content.decode("utf-8").strip().strip('"')
        # s3_uri endpoints use double-encoded s3 uri path parameters
        s3_path_url = quote(quote(s3_path, safe=""))
        schema_res = session.get(f"{api_url}/cellxgene/s3_uri/{s3_path_url}/api/v0.3/schema", allow_redirects=False)
        assert schema_res.status_code == requests.codes.ok
        return schema_res

    s3_uri_response = get_s3_uri()
    return get_schema(s3_uri_response)
