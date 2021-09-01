from backend.corpora.lambdas.api.v1.collection import create_collection
import base64
import json
import os
import time
import unittest
import requests
from requests import HTTPError
from backend.corpora.common.corpora_config import CorporaAuthConfig

API_URL = {
    "prod": "https://api.cellxgene.cziscience.com",
    "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
    "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
    "test": "http://localhost:5000",
}

AUDIENCE = {
    "prod": "cellxgene.cziscience.com/",
    "staging": "cellxgene.staging.single-cell.czi.technology/",
    "test": "cellxgene.dev.single-cell.czi.technology/",
    "dev": "cellxgene.dev.single-cell.czi.technology/",
}


class TestRevisions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.deployment_stage = os.environ["DEPLOYMENT_STAGE"]

    def setUp(self):
        super().setUp()
        self.get_auth_token()
        self.api = API_URL.get(self.deployment_stage)
        self.test_collection_id = "005d611a-14d5-4fbf-846e-571a1f874f70"
        self.test_file_id = "7c93775542b056e048aa474535b8e5c2"
        self.bad_collection_id = "DNE"
        self.bad_file_id = "DNE"

    @classmethod
    def get_auth_token(cls):
        config = CorporaAuthConfig()

        response = requests.post(
            "https://czi-cellxgene-dev.us.auth0.com/oauth/token",
            headers={"content-type": "application/x-www-form-urlencoded"},
            data=dict(
                grant_type="password",
                username=config.test_account_username,
                password=config.test_account_password,
                audience=AUDIENCE.get(cls.deployment_stage),
                scope="openid profile email offline",
                client_id=config.client_id,
                client_secret=config.client_secret,
            ),
        )

        id_token = response.json()["id_token"]
        token = {"id_token": id_token}
        cls.cookie = base64.b64encode(json.dumps(dict(token)).encode("utf-8")).decode()


    def test_version(self):
        res = requests.get(f"{self.api}/dp/v1/deployed_version")
        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)
        self.assertIsNotNone(res.json()["Data Portal"])

    def test_auth(self):
        headers = {"Cookie": f"cxguser={self.cookie}", "Content-Type": "application/json"}
        res = requests.get(f"{self.api}/dp/v1/userinfo", headers=headers)
        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)
        data = json.loads(res.content)
        self.assertEqual(data["email"], "user@example.com")
        self.assertTrue(data["is_authenticated"])

    def test_root_route(self):
        res = requests.get(f"{self.api}/")

        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)

    def test_get_collections(self):
        res = requests.get(f"{self.api}/dp/v1/collections")

        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)
        data = json.loads(res.content)
        for collection in data["collections"]:
            self.assertIsInstance(collection["id"], str)
            self.assertIsInstance(collection["created_at"], float)

    def create_collection(self, headers):
        data = {
            "contact_email": "lisbon@gmail.com",
            "contact_name": "Madrid Sparkle",
            "data_submission_policy_version": "1",
            "description": "Well here are some words",
            "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "protocol.com"}],
            "name": "my2collection",
        }

        res = requests.post(f"{self.api}/dp/v1/collections", data=json.dumps(data), headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        collection_uuid = data["collection_uuid"]
        # Doesn't work since the collection is published. 
        # See https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data-portal/1375
        self.addCleanup(requests.delete, f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
        print(f"Collection id: {collection_uuid}")
        self.assertEqual(res.status_code, requests.codes.created)
        self.assertIn("collection_uuid", data)
        return collection_uuid

    def create_explorer_url(self, dataset_id):
        return f"https://cellxgene.{self.deployment_stage}.single-cell.czi.technology/e/{dataset_id}.cxg/"

    @unittest.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
    def test_revision_flow(self):

        headers = {"Cookie": f"cxguser={self.cookie}", "Content-Type": "application/json"}

        collection_uuid = self.create_collection(headers)

        dataset_1_dropbox_url = "https://www.dropbox.com/s/ib9pth7jr5fvaa8/7MB.h5ad?dl=0"
        dataset_2_dropbox_url = "https://www.dropbox.com/s/ogempho8wixdlvw/local.h5ad?dl=0"

        # Uploads a dataset
        self.upload_and_wait(collection_uuid, dataset_1_dropbox_url)

        # make collection public
        with self.subTest("Test make collection public"):
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}/publish", headers=headers)
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

        dataset_id = requests.get(f"{self.api}/dp/v1/collections/{collection_uuid}").json()['datasets'][0]["id"]
        explorer_url = self.create_explorer_url(dataset_id)

        dataset_meta_payload_before_revision_res = requests.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}")
        dataset_meta_payload_before_revision_res.raise_for_status()
        dataset_meta_payload_before_revision = dataset_meta_payload_before_revision_res.json()

        dataset_schema_before_revision_res = requests.get(f"{self.api}/cellxgene/e/{dataset_id}.cxg/api/v0.2/schema")
        dataset_schema_before_revision_res.raise_for_status()
        dataset_schema_before_revision = dataset_schema_before_revision_res.json()

        with self.subTest("Test updating a dataset in a revision does not effect the published dataset"):
            # Start a revision
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
            self.assertEqual(res.status_code, 201)
            private_dataset_id = res.json()['datasets'][0]['id']

            # Upload a new dataset
            new_dataset_id = self.upload_and_wait(collection_uuid, dataset_2_dropbox_url, existing_dataset_id=private_dataset_id)

            # Check that the published dataset is still the same
            dataset_meta_payload_after_revision = requests.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}").json()
            self.assertDictEqual(dataset_meta_payload_before_revision, dataset_meta_payload_after_revision)
            dataset_schema_after_revision = requests.get(f"{self.api}/cellxgene/e/{dataset_id}.cxg/api/v0.2/schema").json()
            self.assertDictEqual(dataset_schema_before_revision, dataset_schema_after_revision)

        with self.subTest("Publishing a revised dataset replaces the original dataset"):
            # Publish the revision
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}/publish", headers=headers)
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

            dataset_meta_payload = requests.get(f"{self.api}/dp/v1/datasets/meta?url={explorer_url}").json()
            self.assertEqual(dataset_meta_payload["s3_uri"], f"s3://hosted-cellxgene-staging/{new_dataset_id}.cxg/")

            # TODO: add `And the explorer url redirects appropriately`

        with self.subTest("Adding a dataset to a revision does not impact public datasets in that collection"):

            # Get datasets for the collection (before uploading)
            public_datasets_before = requests.get(f"{self.api}/dp/v1/collections/{collection_uuid}").json()["datasets"]

            # Start a revision
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
            self.assertEqual(res.status_code, 201)

            # Upload a new dataset
            another_dataset_id = self.upload_and_wait(collection_uuid, dataset_1_dropbox_url)

            # Get datasets for the collection (after uploading)
            public_datasets_after = requests.get(f"{self.api}/dp/v1/collections/{collection_uuid}").json()["datasets"]
            self.assertCountEqual(public_datasets_before, public_datasets_after)

        with self.subTest("Publishing a revision that contains a new dataset updates the collection page for the data portal (with the new dataset)"):
            # Publish the revision
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}/publish", headers=headers)
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

            # Check if the last updated dataset_id is among the public datasets
            public_datasets = requests.get(f"{self.api}/dp/v1/collections/{collection_uuid}").json()["datasets"]
            self.assertEqual(len(public_datasets), 2)
            ids = [dataset["id"] for dataset in public_datasets]
            self.assertIn(another_dataset_id, ids)

        # deleted_dataset_id = None
        with self.subTest("Deleting a dataset does not effect the published dataset"):
            # Start a revision
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
            self.assertEqual(res.status_code, 201)

            # This only works if you pick the non replaced dataset.
            dataset_to_delete = res.json()['datasets'][1]
            deleted_dataset_id = dataset_to_delete['id']
            original_dataset_id = dataset_to_delete['original_id']
            original_explorer_url = self.create_explorer_url(original_dataset_id)

            # Delete a dataset within the revision
            res = requests.delete(f"{self.api}/dp/v1/datasets/{deleted_dataset_id}", headers=headers)
            self.assertEqual(res.status_code, 202)

            # Check if the dataset is still available
            res = requests.get(f"{self.api}/dp/v1/datasets/meta?url={original_explorer_url}")
            self.assertEqual(res.status_code, 200)

            res = requests.get(f"{self.api}/cellxgene/e/{original_dataset_id}.cxg/api/v0.2/schema")
            self.assertEqual(res.status_code, 200)

        with self.subTest("Publishing a revision that deletes a dataset removes it from the data portal"):
            # Publish the revision
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}/publish", headers=headers)
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

            # Check that the dataset doesn't exist anymore
            res = requests.get(f"{self.api}/dp/v1/datasets/meta?url={self.create_explorer_url(deleted_dataset_id)}")
            self.assertEqual(res.status_code, 404)

            # Deletion happens with eventual consistency, so keeps retrying for 20 seconds
            (final_status_code, desired_status_code) = (None, 404)
            for i in range(20):
                res = requests.get(f"{self.api}/cellxgene/e/{original_dataset_id}.cxg/api/v0.2/schema")
                final_status_code = res.status_code
                if final_status_code == desired_status_code:
                    break
                time.sleep(1)
            self.assertEqual(final_status_code, desired_status_code)

    def upload_and_wait(self, collection_uuid, dropbox_url, existing_dataset_id = None):
        headers = {"Cookie": f"cxguser={self.cookie}", "Content-Type": "application/json"}
        body = {"url": dropbox_url}

        if existing_dataset_id is None:
            res = requests.post(
                f"{self.api}/dp/v1/collections/{collection_uuid}/upload-links", data=json.dumps(body), headers=headers
            )
        else:
            body["id"] = existing_dataset_id
            res = requests.put(
                f"{self.api}/dp/v1/collections/{collection_uuid}/upload-links", data=json.dumps(body), headers=headers
            )

        res.raise_for_status()
        dataset_uuid = json.loads(res.content)["dataset_uuid"]
        self.addCleanup(requests.delete, f"{self.api}/dp/v1/datasets/{dataset_uuid}", headers=headers)

        keep_trying = True
        timer = time.time()
        while keep_trying:
            res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_uuid}/status", headers=headers)
            res.raise_for_status()
            data = json.loads(res.content)
            upload_status = data["upload_status"]
            if upload_status == "UPLOADED":
                conversion_cxg_status = data.get("conversion_cxg_status")
                conversion_loom_status = data.get("conversion_loom_status")
                conversion_rds_status = data.get("conversion_rds_status")
                conversion_anndata_status = data.get("conversion_anndata_status")
                if (
                    conversion_cxg_status
                    == conversion_loom_status
                    == conversion_rds_status
                    == conversion_anndata_status
                    == "CONVERTED"
                ):
                    keep_trying = False
            if time.time() >= timer + 300:
                raise TimeoutError(
                    f"Dataset upload or conversion timed out after 5 min. Check logs for dataset: {dataset_uuid}"
                )
            time.sleep(10)
        return dataset_uuid
