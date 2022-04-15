import json
import os
import unittest
import requests
from tenacity import retry, stop_after_attempt, wait_fixed

from tests.functional.backend.corpora.common import BaseFunctionalTestCase


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
            "description": "Well here are some words",
            "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "protocol.com"}],
            "name": "my2collection",
        }

        res = requests.post(f"{self.api}/dp/v1/collections", data=json.dumps(data), headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        collection_uuid = data["collection_uuid"]

        # Doesn't work since the collection is published. See issue #1375
        self.addCleanup(requests.delete, f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
        self.assertEqual(res.status_code, requests.codes.created)
        self.assertIn("collection_uuid", data)
        return collection_uuid

    def create_explorer_url(self, dataset_id):
        return f"https://cellxgene.{self.deployment_stage}.single-cell.czi.technology/e/{dataset_id}.cxg/"

    @unittest.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
    def test_revision_flow(self):

        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}

        collection_uuid = self.create_collection(headers)

        dataset_1_dropbox_url = "https://www.dropbox.com/s/qiclvn1slmap351/example_valid.h5ad?dl=0"
        dataset_2_dropbox_url = "https://www.dropbox.com/s/qiclvn1slmap351/example_valid.h5ad?dl=0"

        # Uploads a dataset
        self.upload_and_wait(collection_uuid, dataset_1_dropbox_url)

        # make collection public
        with self.subTest("Test make collection public"):
            body = {"data_submission_policy_version": "1.0"}
            res = requests.post(
                f"{self.api}/dp/v1/collections/{collection_uuid}/publish", headers=headers, data=json.dumps(body)
            )
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

        dataset_id = requests.get(f"{self.api}/dp/v1/collections/{collection_uuid}").json()["datasets"][0]["id"]
        explorer_url = self.create_explorer_url(dataset_id)

        meta_payload_before_revision_res = requests.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}")
        meta_payload_before_revision_res.raise_for_status()
        meta_payload_before_revision = meta_payload_before_revision_res.json()

        # Endpoint is eventually consistent
        schema_before_revision = self.get_schema_with_retries(dataset_id).json()

        with self.subTest("Test updating a dataset in a revision does not effect the published dataset"):
            # Start a revision
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
            self.assertEqual(res.status_code, 201)
            revision_uuid = res.json()["id"]
            private_dataset_id = res.json()["datasets"][0]["id"]

            meta_payload_res = requests.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}")
            meta_payload_res.raise_for_status()
            meta_payload = meta_payload_res.json()

            self.assertDictEqual(meta_payload_before_revision, meta_payload)

            # Upload a new dataset
            new_dataset_id = self.upload_and_wait(
                revision_uuid,
                dataset_2_dropbox_url,
                existing_dataset_id=private_dataset_id,
            )

            # Check that the published dataset is still the same
            meta_payload_after_revision = requests.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}").json()
            self.assertDictEqual(meta_payload_before_revision, meta_payload_after_revision)
            schema_after_revision = requests.get(f"{self.api}/cellxgene/e/{dataset_id}.cxg/api/v0.2/schema").json()
            self.assertDictEqual(schema_before_revision, schema_after_revision)

        with self.subTest("Publishing a revised dataset replaces the original dataset"):
            # Publish the revision
            body = {"data_submission_policy_version": "1.0"}
            res = requests.post(
                f"{self.api}/dp/v1/collections/{revision_uuid}/publish", headers=headers, data=json.dumps(body)
            )
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

            dataset_meta_payload = requests.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}").json()
            self.assertEqual(
                dataset_meta_payload["s3_uri"],
                f"s3://hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}/{new_dataset_id}.cxg/",
            )

            # TODO: add `And the explorer url redirects appropriately`

        with self.subTest("Adding a dataset to a revision does not impact public datasets in that collection"):

            # Get datasets for the collection (before uploading)
            public_datasets_before = requests.get(f"{self.api}/dp/v1/collections/{collection_uuid}").json()["datasets"]

            # Start a revision
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
            revision_uuid = res.json()["id"]
            self.assertEqual(res.status_code, 201)

            # Upload a new dataset
            another_dataset_id = self.upload_and_wait(revision_uuid, dataset_1_dropbox_url)

            # Get datasets for the collection (after uploading)
            public_datasets_after = requests.get(f"{self.api}/dp/v1/collections/{collection_uuid}").json()["datasets"]
            self.assertCountEqual(public_datasets_before, public_datasets_after)

        with self.subTest(
            "Publishing a revision that contains a new dataset updates "
            "the collection page for the data portal (with the new dataset)"
        ):
            # Publish the revision
            body = {"data_submission_policy_version": "1.0"}
            res = requests.post(
                f"{self.api}/dp/v1/collections/{revision_uuid}/publish", headers=headers, data=json.dumps(body)
            )
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

            # Check if the last updated dataset_id is among the public datasets
            public_datasets = requests.get(f"{self.api}/dp/v1/collections/{collection_uuid}").json()["datasets"]
            self.assertEqual(len(public_datasets), 2)
            ids = [dataset["id"] for dataset in public_datasets]
            self.assertIn(another_dataset_id, ids)

        with self.subTest("Deleting a dataset does not effect the published dataset"):
            # Start a revision
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
            self.assertEqual(res.status_code, 201)
            revision_uuid = res.json()["id"]

            # This only works if you pick the non replaced dataset.
            dataset_to_delete = res.json()["datasets"][1]
            deleted_dataset_id = dataset_to_delete["id"]
            original_dataset_id = dataset_to_delete["original_id"]
            original_explorer_url = self.create_explorer_url(original_dataset_id)

            # Delete a dataset within the revision
            res = requests.delete(f"{self.api}/dp/v1/datasets/{deleted_dataset_id}", headers=headers)
            self.assertEqual(res.status_code, 202)

            # Check if the dataset is still available
            res = requests.get(f"{self.api}/dp/v1/datasets/meta?url={original_explorer_url}")
            self.assertEqual(res.status_code, 200)

            # Endpoint is eventually consistent
            res = self.get_schema_with_retries(original_dataset_id)
            self.assertEqual(res.status_code, 200)

        with self.subTest("Publishing a revision that deletes a dataset removes it from the data portal"):
            # Publish the revision
            body = {"data_submission_policy_version": "1.0"}
            res = requests.post(
                f"{self.api}/dp/v1/collections/{revision_uuid}/publish", headers=headers, data=json.dumps(body)
            )
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

            # Check that the dataset doesn't exist anymore
            res = requests.get(f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
            res.raise_for_status()
            datasets = [dataset["id"] for dataset in res.json()["datasets"]]
            self.assertEqual(1, len(datasets))
            self.assertNotIn(deleted_dataset_id, datasets)
            self.assertNotIn(original_dataset_id, datasets)

            # Endpoint is eventually consistent. This redirects to the collection page, so the status we want is 302
            res = self.get_schema_with_retries(original_dataset_id, desired_http_status_code=302)
            self.assertEqual(res.status_code, 302)

    @retry(wait=wait_fixed(1), stop=stop_after_attempt(50))
    def get_schema_with_retries(self, dataset_id, desired_http_status_code=requests.codes.ok):
        schema_res = requests.get(f"{self.api}/cellxgene/e/{dataset_id}.cxg/api/v0.2/schema", allow_redirects=False)

        if schema_res.status_code != desired_http_status_code:
            raise UndesiredHttpStatusCodeError

        return schema_res
