import json
import os
import unittest
from urllib.parse import quote

import requests
from tenacity import retry, stop_after_attempt, wait_fixed

from backend.common.constants import DATA_SUBMISSION_POLICY_VERSION
from tests.functional.backend.constants import ATAC_SEQ_MANIFEST, DATASET_URI
from tests.functional.backend.skip_reason import skip_creation_on_prod
from tests.functional.backend.utils import assertStatusCode, create_explorer_url, create_test_collection


def _setup_published_collection_with_dataset(
    headers, request, session, api_url, collection_data, upload_func, upload_data
):
    """Helper to create a collection, upload a dataset, and publish it."""
    # Create a test collection
    collection_id = create_test_collection(headers, request, session, api_url, collection_data)

    # Upload dataset using provided function and data
    upload_func(collection_id, upload_data)

    # Make collection public
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(f"{api_url}/dp/v1/collections/{collection_id}/publish", headers=headers, data=json.dumps(body))
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # Get canonical collection ID post-publish
    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    data = json.loads(res.content)
    canonical_collection_id = data["id"]

    return canonical_collection_id, collection_id


def _get_dataset_info_and_baselines(session, api_url, canonical_collection_id):
    """Helper to get dataset info and baseline metadata/schema for comparison."""
    dataset_response = session.get(f"{api_url}/dp/v1/collections/{canonical_collection_id}").json()["datasets"][0]
    dataset_id = dataset_response["id"]
    explorer_url = dataset_response["dataset_deployments"][0]["url"]

    # Get baseline metadata and schema
    if os.environ["DEPLOYMENT_STAGE"] == "rdev":
        explorer_url = explorer_url.replace("rdev.", "")
        meta_payload_before_revision = {}
        schema_before_revision = {}
    else:
        meta_payload_before_revision_res = session.get(f"{api_url}/dp/v1/datasets/meta?url={explorer_url}")
        meta_payload_before_revision_res.raise_for_status()
        meta_payload_before_revision = meta_payload_before_revision_res.json()
        schema_before_revision = get_schema_with_retries(dataset_id, api_url, session).json()

    return dataset_id, explorer_url, meta_payload_before_revision, schema_before_revision


def _test_dataset_replacement_flow(
    session,
    api_url,
    headers,
    canonical_collection_id,
    explorer_url,
    meta_payload_before_revision,
    schema_before_revision,
    upload_func,
    upload_data,
):
    """Helper to test dataset replacement via revision."""
    # Start a revision
    res = session.post(f"{api_url}/dp/v1/collections/{canonical_collection_id}", headers=headers)
    assertStatusCode(201, res)
    data = json.loads(res.content)
    revision_id = data["id"]
    private_dataset_id = res.json()["datasets"][0]["id"]

    # Verify published dataset is unchanged after revision creation
    if os.environ["DEPLOYMENT_STAGE"] != "rdev":
        meta_payload_res = session.get(f"{api_url}/dp/v1/datasets/meta?url={explorer_url}")
        meta_payload_res.raise_for_status()
        meta_payload = meta_payload_res.json()
        assert meta_payload_before_revision == meta_payload

    # Upload replacement dataset
    upload_func(
        revision_id,
        upload_data,
        existing_dataset_id=private_dataset_id,
    )

    # Verify published dataset is still unchanged after upload
    if os.environ["DEPLOYMENT_STAGE"] != "rdev":
        meta_payload_after_revision = session.get(f"{api_url}/dp/v1/datasets/meta?url={explorer_url}").json()
        assert meta_payload_before_revision == meta_payload_after_revision
        # Get the original dataset ID to check schema (need to extract from canonical collection)
        original_dataset_id = session.get(f"{api_url}/dp/v1/collections/{canonical_collection_id}").json()["datasets"][
            0
        ]["id"]
        schema_after_revision = get_schema_with_retries(original_dataset_id, api_url, session).json()
        assert schema_before_revision == schema_after_revision

    # Publish the revision
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(f"{api_url}/dp/v1/collections/{revision_id}/publish", headers=headers, data=json.dumps(body))
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # Verify the dataset was replaced
    if os.environ["DEPLOYMENT_STAGE"] != "rdev":
        dataset_meta_payload = session.get(f"{api_url}/dp/v1/datasets/meta?url={explorer_url}").json()
        assert dataset_meta_payload["s3_uri"].startswith(f"s3://hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}/")
        assert dataset_meta_payload["s3_uri"].endswith(".cxg/")
        assert (
            dataset_meta_payload["dataset_id"] in dataset_meta_payload["s3_uri"]
        ), "The S3_URI should contain the revised dataset id."


def _test_dataset_addition_flow(session, api_url, headers, canonical_collection_id, upload_func, upload_data):
    """Helper to test adding new datasets via revision."""
    # Start a new revision
    res = session.post(f"{api_url}/dp/v1/collections/{canonical_collection_id}", headers=headers)
    assertStatusCode(201, res)
    revision_id = res.json()["id"]

    # Get datasets before uploading
    public_datasets_before = session.get(f"{api_url}/dp/v1/collections/{canonical_collection_id}").json()["datasets"]

    # Upload a new dataset
    another_dataset_result = upload_func(revision_id, upload_data)
    another_dataset_id = another_dataset_result["dataset_id"]

    # Verify adding dataset to revision doesn't affect public datasets
    public_datasets_after = session.get(f"{api_url}/dp/v1/collections/{canonical_collection_id}").json()["datasets"]
    unittest.TestCase().assertCountEqual(public_datasets_before, public_datasets_after)

    # Publish the revision
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(f"{api_url}/dp/v1/collections/{revision_id}/publish", headers=headers, data=json.dumps(body))
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # Verify the new dataset is now public
    public_datasets = session.get(f"{api_url}/dp/v1/collections/{canonical_collection_id}").json()["datasets"]
    assert len(public_datasets) == 2
    ids = [dataset["id"] for dataset in public_datasets]
    assert another_dataset_id in ids


def _test_dataset_deletion_flow(
    session, api_url, headers, canonical_collection_id, collection_id, curation_api_access_token, deployment_stage
):
    """Helper to test deleting datasets via revision."""
    # Start a revision
    res = session.post(f"{api_url}/dp/v1/collections/{canonical_collection_id}", headers=headers)
    assertStatusCode(201, res)
    revision_id = res.json()["id"]

    # Pick the second dataset to delete (the one added in the previous step)
    dataset_to_delete = res.json()["datasets"][1]
    revision_deleted_dataset_id = dataset_to_delete["id"]
    published_explorer_url = create_explorer_url(revision_deleted_dataset_id, deployment_stage)

    # Delete the dataset from the revision
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

    # Verify deleting in revision doesn't affect published dataset
    res = session.get(f"{api_url}/dp/v1/datasets/meta?url={published_explorer_url}")
    assertStatusCode(200, res)

    res = get_schema_with_retries(revision_deleted_dataset_id, api_url, session)
    assertStatusCode(200, res)

    # Publish the revision
    body = {"data_submission_policy_version": DATA_SUBMISSION_POLICY_VERSION}
    res = session.post(f"{api_url}/dp/v1/collections/{revision_id}/publish", headers=headers, data=json.dumps(body))
    res.raise_for_status()
    assertStatusCode(requests.codes.accepted, res)

    # Verify the dataset is now deleted from the collection
    res = session.get(f"{api_url}/dp/v1/collections/{collection_id}", headers=headers)
    res.raise_for_status()
    datasets = [dataset["id"] for dataset in res.json()["datasets"]]
    assert len(datasets) == 1
    assert revision_deleted_dataset_id not in datasets


@skip_creation_on_prod
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
    """Test revision flow using dataset uploaded via URL."""
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}

    # Setup: create collection, upload dataset, and publish
    canonical_collection_id, collection_id = _setup_published_collection_with_dataset(
        headers, request, session, api_url, collection_data, upload_dataset, DATASET_URI
    )

    # Get dataset info and baseline metadata/schema
    dataset_id, explorer_url, meta_payload_before_revision, schema_before_revision = _get_dataset_info_and_baselines(
        session, api_url, canonical_collection_id
    )

    # Test dataset replacement flow
    _test_dataset_replacement_flow(
        session,
        api_url,
        headers,
        canonical_collection_id,
        explorer_url,
        meta_payload_before_revision,
        schema_before_revision,
        upload_dataset,
        DATASET_URI,
    )

    # Test dataset addition flow
    _test_dataset_addition_flow(session, api_url, headers, canonical_collection_id, upload_dataset, DATASET_URI)

    # Test dataset deletion flow
    _test_dataset_deletion_flow(
        session, api_url, headers, canonical_collection_id, collection_id, curation_api_access_token, deployment_stage
    )


@skip_creation_on_prod
def test_atac_revision_flow(
    curator_cookie,
    session,
    api_url,
    upload_manifest,
    curation_api_access_token,
    deployment_stage,
    request,
    collection_data,
):
    """Test revision flow using ATAC dataset uploaded via manifest."""
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}

    # Setup: create collection, upload ATAC dataset, and publish
    canonical_collection_id, collection_id = _setup_published_collection_with_dataset(
        headers, request, session, api_url, collection_data, upload_manifest, ATAC_SEQ_MANIFEST
    )

    # Get dataset info and baseline metadata/schema
    dataset_id, explorer_url, meta_payload_before_revision, schema_before_revision = _get_dataset_info_and_baselines(
        session, api_url, canonical_collection_id
    )

    # Test ATAC dataset replacement flow
    _test_dataset_replacement_flow(
        session,
        api_url,
        headers,
        canonical_collection_id,
        explorer_url,
        meta_payload_before_revision,
        schema_before_revision,
        upload_manifest,
        ATAC_SEQ_MANIFEST,
    )

    # Test ATAC dataset addition flow
    _test_dataset_addition_flow(session, api_url, headers, canonical_collection_id, upload_manifest, ATAC_SEQ_MANIFEST)

    # Test ATAC dataset deletion flow
    _test_dataset_deletion_flow(
        session, api_url, headers, canonical_collection_id, collection_id, curation_api_access_token, deployment_stage
    )


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
