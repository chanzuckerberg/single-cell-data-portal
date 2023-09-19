import json
import os
import unittest
from urllib.parse import quote

import requests
from tenacity import retry, stop_after_attempt, wait_fixed

from tests.functional.backend.common import BaseFunctionalTestCase


class UndesiredHttpStatusCodeError(Exception):
    pass


class TestRevisions(BaseFunctionalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

    def create_collection(self, headers):
        data = {
            "contact_email": "lisbon@gmail.com",
            "contact_name": "Madrid Sparkle",
            "curator_name": "John Smith",
            "description": "Well here are some words",
            "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "http://protocol.com"}],
            "name": "my2collection",
        }

        res = self.session.post(f"{self.api}/dp/v1/collections", data=json.dumps(data), headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        collection_id = data["collection_id"]

        # Doesn't work since the collection is published. See issue #1375
        # Should work now via cxg-admin role thru Curation API
        curation_api_headers = {"Authorization": f"Bearer {self.curation_api_access_token}"}
        self.addCleanup(
            self.session.delete,
            f"{self.api}/curation/v1/collections/{collection_id}?delete_published=true",
            headers=curation_api_headers,
        )
        self.assertStatusCode(requests.codes.created, res)
        self.assertIn("collection_id", data)
        return collection_id

    def create_explorer_url(self, dataset_id):
        return f"https://cellxgene.{self.deployment_stage}.single-cell.czi.technology/e/{dataset_id}.cxg/"

    @unittest.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
    def test_revision_flow(self):

        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}

        collection_id = self.create_collection(headers)

        dataset_1_dropbox_url = "https://www.dropbox.com/s/m1ur46nleit8l3w/3_0_0_valid.h5ad?dl=0"
        dataset_2_dropbox_url = "https://www.dropbox.com/s/m1ur46nleit8l3w/3_0_0_valid.h5ad?dl=0"

        # Uploads a dataset
        self.upload_and_wait(collection_id, dataset_1_dropbox_url)

        # make collection public
        with self.subTest("Test make collection public"):
            body = {"data_submission_policy_version": "1.0"}
            res = self.session.post(
                f"{self.api}/dp/v1/collections/{collection_id}/publish", headers=headers, data=json.dumps(body)
            )
            res.raise_for_status()
            self.assertStatusCode(requests.codes.accepted, res)

        # get canonical collection id, post-publish
        res = self.session.get(f"{self.api}/dp/v1/collections/{collection_id}", headers=headers)
        data = json.loads(res.content)
        canonical_collection_id = data["id"]

        dataset_response = self.session.get(f"{self.api}/dp/v1/collections/{canonical_collection_id}").json()[
            "datasets"
        ][0]
        dataset_id = dataset_response["id"]
        explorer_url = dataset_response["dataset_deployments"][0]["url"]

        meta_payload_before_revision_res = self.session.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}")
        meta_payload_before_revision_res.raise_for_status()
        meta_payload_before_revision = meta_payload_before_revision_res.json()

        # Endpoint is eventually consistent
        schema_before_revision = self.get_schema_with_retries(dataset_id).json()

        # Start a revision
        res = self.session.post(f"{self.api}/dp/v1/collections/{canonical_collection_id}", headers=headers)
        self.assertStatusCode(201, res)
        data = json.loads(res.content)
        revision_id = data["id"]

        with self.subTest("Test updating a dataset in a revision does not effect the published dataset"):
            private_dataset_id = res.json()["datasets"][0]["id"]

            meta_payload_res = self.session.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}")
            meta_payload_res.raise_for_status()
            meta_payload = meta_payload_res.json()

            self.assertDictEqual(meta_payload_before_revision, meta_payload)

            # Upload a new dataset
            self.upload_and_wait(
                revision_id,
                dataset_2_dropbox_url,
                existing_dataset_id=private_dataset_id,
            )

            # Check that the published dataset is still the same
            meta_payload_after_revision = self.session.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}").json()
            self.assertDictEqual(meta_payload_before_revision, meta_payload_after_revision)
            schema_after_revision = self.get_schema_with_retries(dataset_id).json()
            self.assertDictEqual(schema_before_revision, schema_after_revision)

        with self.subTest("Publishing a revised dataset replaces the original dataset"):
            # Publish the revision
            body = {"data_submission_policy_version": "1.0"}
            res = self.session.post(
                f"{self.api}/dp/v1/collections/{revision_id}/publish", headers=headers, data=json.dumps(body)
            )
            res.raise_for_status()
            self.assertStatusCode(requests.codes.accepted, res)

            dataset_meta_payload = self.session.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}").json()
            self.assertTrue(
                dataset_meta_payload["s3_uri"].startswith(f"s3://hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}/")
            )
            self.assertTrue(dataset_meta_payload["s3_uri"].endswith(".cxg/"))
            self.assertIn(
                dataset_meta_payload["dataset_id"],
                dataset_meta_payload["s3_uri"],
                "The id of the S3_URI should be the revised dataset id.",
            )

            # TODO: add `And the explorer url redirects appropriately`

        # Start a new revision
        res = self.session.post(f"{self.api}/dp/v1/collections/{canonical_collection_id}", headers=headers)
        self.assertStatusCode(201, res)
        revision_id = res.json()["id"]

        # Get datasets for the collection (before uploading)
        public_datasets_before = self.session.get(f"{self.api}/dp/v1/collections/{canonical_collection_id}").json()[
            "datasets"
        ]

        # Upload a new dataset
        another_dataset_id = self.upload_and_wait(revision_id, dataset_1_dropbox_url)

        with self.subTest("Adding a dataset to a revision does not impact public datasets in that collection"):
            # Get datasets for the collection (after uploading)
            public_datasets_after = self.session.get(f"{self.api}/dp/v1/collections/{canonical_collection_id}").json()[
                "datasets"
            ]
            self.assertCountEqual(public_datasets_before, public_datasets_after)

        # Publish the revision
        body = {"data_submission_policy_version": "1.0"}
        res = self.session.post(
            f"{self.api}/dp/v1/collections/{revision_id}/publish", headers=headers, data=json.dumps(body)
        )
        res.raise_for_status()
        self.assertStatusCode(requests.codes.accepted, res)

        with self.subTest(
            "Publishing a revision that contains a new dataset updates "
            "the collection page for the data portal (with the new dataset)"
        ):
            # Check if the last updated dataset_id is among the public datasets
            public_datasets = self.session.get(f"{self.api}/dp/v1/collections/{canonical_collection_id}").json()[
                "datasets"
            ]
            self.assertEqual(len(public_datasets), 2)
            ids = [dataset["id"] for dataset in public_datasets]
            self.assertIn(another_dataset_id, ids)

        # Start a revision
        res = self.session.post(f"{self.api}/dp/v1/collections/{canonical_collection_id}", headers=headers)
        self.assertStatusCode(201, res)
        revision_id = res.json()["id"]

        # This only works if you pick the non replaced dataset.
        dataset_to_delete = res.json()["datasets"][1]
        revision_deleted_dataset_id = dataset_to_delete["id"]
        published_explorer_url = self.create_explorer_url(revision_deleted_dataset_id)

        # Delete (tombstone) a dataset (using admin privileges) within the revision
        revision_datasets = self.session.get(f"{self.api}/curation/v1/collections/{revision_id}").json()["datasets"]
        dataset_id_to_delete = None
        for dataset in revision_datasets:
            if dataset["dataset_version_id"] == revision_deleted_dataset_id:
                dataset_id_to_delete = dataset["dataset_id"]

        curation_api_headers = {"Authorization": f"Bearer {self.curation_api_access_token}"}
        res = self.session.delete(
            f"{self.api}/curation/v1/collections/{revision_id}/datasets/{dataset_id_to_delete}?delete_published=true",
            headers=curation_api_headers,
        )
        self.assertStatusCode(202, res)

        with self.subTest("Deleting a dataset does not effect the published dataset"):
            # Check if the dataset is still available
            res = self.session.get(f"{self.api}/dp/v1/datasets/meta?url={published_explorer_url}")
            self.assertStatusCode(200, res)

            # Endpoint is eventually consistent
            res = self.get_schema_with_retries(revision_deleted_dataset_id)
            self.assertStatusCode(200, res)

        with self.subTest("Publishing a revision that deletes a dataset removes it from the data portal"):
            # Publish the revision
            body = {"data_submission_policy_version": "1.0"}
            res = self.session.post(
                f"{self.api}/dp/v1/collections/{revision_id}/publish", headers=headers, data=json.dumps(body)
            )
            res.raise_for_status()
            self.assertStatusCode(requests.codes.accepted, res)

            # Check that the dataset doesn't exist anymore
            res = self.session.get(f"{self.api}/dp/v1/collections/{collection_id}", headers=headers)
            res.raise_for_status()
            datasets = [dataset["id"] for dataset in res.json()["datasets"]]
            self.assertEqual(1, len(datasets))
            self.assertNotIn(revision_deleted_dataset_id, datasets)

    def get_schema_with_retries(self, dataset_id, desired_http_status_code=requests.codes.ok):
        @retry(wait=wait_fixed(1), stop=stop_after_attempt(50))
        def get_s3_uri():
            s3_uri_res = self.session.get(
                f"{self.api}/cellxgene/e/{dataset_id}.cxg/api/v0.3/s3_uri", allow_redirects=False
            )
            if s3_uri_res.status_code != desired_http_status_code:
                raise UndesiredHttpStatusCodeError
            return s3_uri_res

        @retry(wait=wait_fixed(1), stop=stop_after_attempt(50))
        def get_schema(s3_uri_response_object):
            # parse s3_uri_response_object content
            s3_path = s3_uri_response_object.content.decode("utf-8").strip().strip('"')
            # s3_uri endpoints use double-encoded s3 uri path parameters
            s3_path_url = quote(quote(s3_path, safe=""))
            schema_res = self.session.get(
                f"{self.api}/cellxgene/s3_uri/{s3_path_url}/api/v0.3/schema", allow_redirects=False
            )
            if schema_res.status_code != requests.codes.ok:
                raise UndesiredHttpStatusCodeError
            return schema_res

        s3_uri_response = get_s3_uri()
        return get_schema(s3_uri_response)
