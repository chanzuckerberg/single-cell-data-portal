import json
import os
import time
import unittest

import requests
from requests import HTTPError

from tests.functional.backend.common import TEST_DATASET_URI, BaseFunctionalTestCase


class TestApi(BaseFunctionalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

    @unittest.skip("Skipping this test until the deployed_version feature is fixed")
    def test_version(self):
        res = self.session.get(f"{self.api}/dp/v1/deployed_version")
        res.raise_for_status()
        self.assertStatusCode(requests.codes.ok, res)
        self.assertTrue(len(res.json()["Data Portal"]) > 0)

    def test_auth(self):
        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}
        res = self.session.get(f"{self.api}/dp/v1/userinfo", headers=headers)
        res.raise_for_status()
        self.assertStatusCode(requests.codes.ok, res)
        data = json.loads(res.content)
        self.assertEqual(data["email"], "functest@example.com")

    def test_root_route(self):
        res = self.session.get(f"{self.api}/")

        res.raise_for_status()
        self.assertStatusCode(requests.codes.ok, res)

    def test_get_collections(self):
        res = self.session.get(f"{self.api}/dp/v1/collections")

        res.raise_for_status()
        self.assertStatusCode(requests.codes.ok, res)
        data = json.loads(res.content)
        for collection in data["collections"]:
            self.assertIsInstance(collection["id"], str)
            self.assertIsInstance(collection["created_at"], float)

    @unittest.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
    def test_collection_flow(self):
        # create collection
        data = {
            "contact_email": "lisbon@gmail.com",
            "contact_name": "Madrid Sparkle",
            "curator_name": "John Smith",
            "description": "Well here are some words",
            "links": [
                {"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "https://protocol.com"}
            ],
            "name": "my2collection",
        }

        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}
        res = self.session.post(f"{self.api}/dp/v1/collections", data=json.dumps(data), headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        collection_id = data["collection_id"]
        self.assertStatusCode(requests.codes.created, res)
        self.assertIn("collection_id", data)

        with self.subTest("Test created collection is private"):
            res = self.session.get(f"{self.api}/dp/v1/collections", headers=headers)
            data = json.loads(res.content)
            private_collection_ids = []
            for collection in data["collections"]:
                if collection["visibility"] == "PRIVATE":
                    private_collection_ids.append(collection["id"])
            self.assertIn(collection_id, private_collection_ids)

        with self.subTest("Test update collection info"):
            updated_data = {
                "contact_email": "person@random.com",
                "contact_name": "Doctor Who",
                "description": "These are different words",
                "links": [
                    {"link_name": "The Source", "link_type": "DATA_SOURCE", "link_url": "https://datasource.com"}
                ],
                "name": "lots of cells",
            }
            res = self.session.put(
                f"{self.api}/dp/v1/collections/{collection_id}", data=json.dumps(updated_data), headers=headers
            )
            res.raise_for_status()
            data = json.loads(res.content)
            data.pop("access_type")
            for key in updated_data:
                self.assertEqual(updated_data[key], data[key])

        self.upload_and_wait(collection_id, TEST_DATASET_URI)

        # make collection public
        with self.subTest("Test make collection public"):
            body = {"data_submission_policy_version": "1.0"}
            res = self.session.post(
                f"{self.api}/dp/v1/collections/{collection_id}/publish", headers=headers, data=json.dumps(body)
            )
            res.raise_for_status()
            self.assertStatusCode(requests.codes.accepted, res)

            # get canonical collection_id
            res = self.session.get(f"{self.api}/dp/v1/collections/{collection_id}", headers=headers)
            data = json.loads(res.content)
            canonical_collection_id = data["id"]

            # check collection returns as public
            res = self.session.get(f"{self.api}/dp/v1/collections", headers=headers)
            data = json.loads(res.content)
            public_collection_ids = []
            for collection in data["collections"]:
                if collection["visibility"] == "PUBLIC":
                    public_collection_ids.append(collection["id"])

            self.assertIn(canonical_collection_id, public_collection_ids)

        with self.subTest("Test everyone can retrieve a public collection"):
            no_auth_headers = {"Content-Type": "application/json"}
            res = self.session.get(f"{self.api}/dp/v1/collections", headers=no_auth_headers)
            data = json.loads(res.content)
            collection_ids = [x["id"] for x in data["collections"]]
            self.assertIn(canonical_collection_id, collection_ids)

        with self.subTest("Test a public collection cannot be tombstoned"):
            res = self.session.delete(f"{self.api}/dp/v1/collections/{canonical_collection_id}", headers=headers)
            self.assertStatusCode(requests.codes.method_not_allowed, res)

            res = self.session.get(f"{self.api}/dp/v1/collections/{collection_id}", headers=headers)
            self.assertStatusCode(requests.codes.ok, res)

    @unittest.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
    def test_delete_private_collection(self):
        # create collection
        data = {
            "contact_email": "lisbon@gmail.com",
            "contact_name": "Madrid Sparkle",
            "curator_name": "John Smith",
            "description": "Well here are some words",
            "links": [
                {"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "https://protocol.com"}
            ],
            "name": "my2collection",
        }

        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}
        res = self.session.post(f"{self.api}/dp/v1/collections", data=json.dumps(data), headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        collection_id = data["collection_id"]
        self.addCleanup(self.session.delete, f"{self.api}/dp/v1/collections/{collection_id}", headers=headers)
        self.assertStatusCode(requests.codes.created, res)
        self.assertIn("collection_id", data)

        # check created collection returns as private
        res = self.session.get(f"{self.api}/dp/v1/collections", headers=headers)
        data = json.loads(res.content)
        private_collection_ids = []
        for collection in data["collections"]:
            if collection["visibility"] == "PRIVATE":
                private_collection_ids.append(collection["id"])
        self.assertIn(collection_id, private_collection_ids)

        # delete collection
        res = self.session.delete(f"{self.api}/dp/v1/collections/{collection_id}?visibility=PRIVATE", headers=headers)
        res.raise_for_status()
        self.assertStatusCode(requests.codes.no_content, res)

        # check collection gone
        no_auth_headers = {"Content-Type": "application/json"}
        res = self.session.get(f"{self.api}/dp/v1/collections?visibility=PRIVATE", headers=no_auth_headers)
        data = json.loads(res.content)
        collection_ids = [x["id"] for x in data["collections"]]
        self.assertNotIn(collection_id, collection_ids)

    @unittest.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
    def test_dataset_upload_flow(self):
        body = {
            "contact_email": "lisbon@gmail.com",
            "contact_name": "Madrid Sparkle",
            "curator_name": "John Smith",
            "description": "Well here are some words",
            "links": [
                {"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "https://protocol.com"}
            ],
            "name": "my2collection",
        }

        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}
        res = self.session.post(f"{self.api}/dp/v1/collections", data=json.dumps(body), headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        collection_id = data["collection_id"]
        self.addCleanup(self.session.delete, f"{self.api}/dp/v1/collections/{collection_id}", headers=headers)
        self.assertStatusCode(requests.codes.created, res)
        self.assertIn("collection_id", data)

        body = {"url": TEST_DATASET_URI}

        res = self.session.post(
            f"{self.api}/dp/v1/collections/{collection_id}/upload-links", data=json.dumps(body), headers=headers
        )
        res.raise_for_status()
        dataset_id = json.loads(res.content)["dataset_id"]
        self.addCleanup(self.session.delete, f"{self.api}/dp/v1/datasets/{dataset_id}", headers=headers)

        self.assertStatusCode(requests.codes.accepted, res)

        res = self.session.get(f"{self.api}/dp/v1/datasets/{dataset_id}/status", headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        self.assertStatusCode(requests.codes.ok, res)
        self.assertEqual(data["upload_status"], "WAITING")

        with self.subTest("Test dataset conversion"):
            keep_trying = True
            expected_upload_statuses = ["WAITING", "UPLOADING", "UPLOADED"]
            # conversion statuses can be `None` when/if we hit the status endpoint too early after an upload
            expected_conversion_statuses = ["CONVERTING", "CONVERTED", "FAILED", "UPLOADING", "UPLOADED", "NA", None]
            timer = time.time()
            while keep_trying:
                data = None
                res = self.session.get(f"{self.api}/dp/v1/datasets/{dataset_id}/status", headers=headers)
                res.raise_for_status()
                data = json.loads(res.content)
                upload_status = data["upload_status"]
                if upload_status:
                    self.assertIn(upload_status, expected_upload_statuses)

                # conversion statuses only returned once uploaded
                if upload_status == "UPLOADED":
                    cxg_status = data.get("cxg_status")
                    rds_status = data.get("rds_status")
                    h5ad_status = data.get("h5ad_status")
                    self.assertIn(data.get("cxg_status"), expected_conversion_statuses)
                    if cxg_status == "FAILED":
                        self.fail(f"CXG CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
                    if rds_status == "FAILED":
                        self.fail(f"RDS CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
                    if h5ad_status == "FAILED":
                        self.fail(f"Anndata CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_id}")
                    if cxg_status == rds_status == h5ad_status == "UPLOADED":
                        keep_trying = False
                if time.time() >= timer + 600:
                    raise TimeoutError(
                        f"Dataset upload or conversion timed out after 10 min. Check logs for dataset: {dataset_id}"
                    )
                time.sleep(10)

        with self.subTest("test non owner cant retrieve status"):
            no_auth_headers = {"Content-Type": "application/json"}
            res = self.session.get(f"{self.api}/dp/v1/datasets/{dataset_id}/status", headers=no_auth_headers)
            with self.assertRaises(HTTPError):
                res.raise_for_status()

        with self.subTest("Test dataset deletion"):
            res = self.session.delete(f"{self.api}/dp/v1/datasets/{dataset_id}", headers=headers)
            res.raise_for_status()
            self.assertStatusCode(requests.codes.accepted, res)

            # Check that the dataset is gone from collection version
            res = self.session.get(f"{self.api}/dp/v1/collections/{collection_id}", headers=headers)
            data = json.loads(res.content)
            datasets = data["datasets"]
            dataset_ids = [dataset.get("id") for dataset in datasets]
            self.assertNotIn(dataset_id, dataset_ids)
