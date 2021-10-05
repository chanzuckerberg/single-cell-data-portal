import json
import os
import time
import unittest
import requests
from requests import HTTPError

from tests.functional.backend.corpora.common import BaseFunctionalTestCase


class TestApi(BaseFunctionalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

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

    @unittest.skipIf(os.environ["DEPLOYMENT_STAGE"] == "prod", "Do not make test collections public in prod")
    def test_collection_flow(self):
        # create collection
        data = {
            "contact_email": "lisbon@gmail.com",
            "contact_name": "Madrid Sparkle",
            "data_submission_policy_version": "1",
            "description": "Well here are some words",
            "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "protocol.com"}],
            "name": "my2collection",
        }

        headers = {"Cookie": f"cxguser={self.cookie}", "Content-Type": "application/json"}
        res = requests.post(f"{self.api}/dp/v1/collections", data=json.dumps(data), headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        collection_uuid = data["collection_uuid"]
        self.assertEqual(res.status_code, requests.codes.created)
        self.assertIn("collection_uuid", data)

        with self.subTest("Test created collection is private"):
            res = requests.get(f"{self.api}/dp/v1/collections", headers=headers)
            data = json.loads(res.content)
            private_collection_uuids = []
            for collection in data["collections"]:
                if collection["visibility"] == "PRIVATE":
                    private_collection_uuids.append(collection["id"])
            self.assertIn(collection_uuid, private_collection_uuids)

        with self.subTest("Test update collection info"):
            updated_data = {
                "contact_email": "person@random.com",
                "contact_name": "Doctor Who",
                "description": "These are different words",
                "links": [{"link_name": "The Source", "link_type": "DATA_SOURCE", "link_url": "datasource.com"}],
                "name": "lots of cells",
            }
            res = requests.put(
                f"{self.api}/dp/v1/collections/{collection_uuid}", data=json.dumps(updated_data), headers=headers
            )
            res.raise_for_status()
            data = json.loads(res.content)
            data.pop("access_type")
            for key in updated_data.keys():
                self.assertEqual(updated_data[key], data[key])

        self.upload_and_wait(collection_uuid, "https://www.dropbox.com/s/qiclvn1slmap351/example_valid.h5ad?dl=0")

        # make collection public
        with self.subTest("Test make collection public"):
            res = requests.post(f"{self.api}/dp/v1/collections/{collection_uuid}/publish", headers=headers)
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

            # check  collection returns as public
            res = requests.get(f"{self.api}/dp/v1/collections", headers=headers)
            data = json.loads(res.content)
            public_collection_uuids = []
            for collection in data["collections"]:
                if collection["visibility"] == "PUBLIC":
                    public_collection_uuids.append(collection["id"])

            self.assertIn(collection_uuid, public_collection_uuids)

        with self.subTest("Test everyone can retrieve a public collection"):
            no_auth_headers = {"Content-Type": "application/json"}
            res = requests.get(f"{self.api}/dp/v1/collections", headers=no_auth_headers)
            data = json.loads(res.content)
            collection_uuids = [x["id"] for x in data["collections"]]
            self.assertIn(collection_uuid, collection_uuids)

        # cannot delete public collection
        with self.subTest("Test a public collection can not be deleted"):
            res = requests.delete(f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
            self.assertEqual(res.status_code, requests.codes.not_allowed)

    def test_delete_private_collection(self):
        # create collection
        data = {
            "contact_email": "lisbon@gmail.com",
            "contact_name": "Madrid Sparkle",
            "data_submission_policy_version": "1",
            "description": "Well here are some words",
            "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "protocol.com"}],
            "name": "my2collection",
        }

        headers = {"Cookie": f"cxguser={self.cookie}", "Content-Type": "application/json"}
        res = requests.post(f"{self.api}/dp/v1/collections", data=json.dumps(data), headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        collection_uuid = data["collection_uuid"]
        self.addCleanup(requests.delete, f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
        self.assertEqual(res.status_code, requests.codes.created)
        self.assertIn("collection_uuid", data)

        # check created collection returns as private
        res = requests.get(f"{self.api}/dp/v1/collections", headers=headers)
        data = json.loads(res.content)
        private_collection_uuids = []
        for collection in data["collections"]:
            if collection["visibility"] == "PRIVATE":
                private_collection_uuids.append(collection["id"])
        self.assertIn(collection_uuid, private_collection_uuids)

        # delete collection
        res = requests.delete(f"{self.api}/dp/v1/collections/{collection_uuid}?visibility=PRIVATE", headers=headers)
        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.no_content)

        # check collection gone
        no_auth_headers = {"Content-Type": "application/json"}
        res = requests.get(f"{self.api}/dp/v1/collections?visibility=PRIVATE", headers=no_auth_headers)
        data = json.loads(res.content)
        collection_uuids = [x["id"] for x in data["collections"]]
        self.assertNotIn(collection_uuid, collection_uuids)

    def test_dataset_upload_flow(self):
        body = {
            "contact_email": "lisbon@gmail.com",
            "contact_name": "Madrid Sparkle",
            "data_submission_policy_version": "1",
            "description": "Well here are some words",
            "links": [{"link_name": "a link to somewhere", "link_type": "PROTOCOL", "link_url": "protocol.com"}],
            "name": "my2collection",
        }

        headers = {"Cookie": f"cxguser={self.cookie}", "Content-Type": "application/json"}
        res = requests.post(f"{self.api}/dp/v1/collections", data=json.dumps(body), headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        collection_uuid = data["collection_uuid"]
        self.addCleanup(requests.delete, f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
        self.assertEqual(res.status_code, requests.codes.created)
        self.assertIn("collection_uuid", data)

        body = {"url": "https://www.dropbox.com/s/qiclvn1slmap351/example_valid.h5ad?dl=0"}

        res = requests.post(
            f"{self.api}/dp/v1/collections/{collection_uuid}/upload-links", data=json.dumps(body), headers=headers
        )
        res.raise_for_status()
        dataset_uuid = json.loads(res.content)["dataset_uuid"]
        self.addCleanup(requests.delete, f"{self.api}/dp/v1/datasets/{dataset_uuid}", headers=headers)

        self.assertEqual(res.status_code, requests.codes.accepted)

        res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_uuid}/status", headers=headers)
        res.raise_for_status()
        data = json.loads(res.content)
        self.assertEqual(res.status_code, requests.codes.ok)
        self.assertEqual(data["upload_status"], "WAITING")

        with self.subTest("Test dataset conversion"):
            keep_trying = True
            expected_upload_statuses = ["WAITING", "UPLOADING", "UPLOADED"]
            # conversion statuses can be `None` when/if we hit the status endpoint too early after an upload
            expected_conversion_statuses = ["CONVERTING", "CONVERTED", "FAILED", None]
            timer = time.time()
            while keep_trying:
                data = None
                res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_uuid}/status", headers=headers)
                res.raise_for_status()
                data = json.loads(res.content)
                upload_status = data["upload_status"]
                if upload_status:
                    self.assertIn(upload_status, expected_upload_statuses)

                # conversion statuses only returned once uploaded
                if upload_status == "UPLOADED":
                    conversion_cxg_status = data.get("conversion_cxg_status")
                    conversion_loom_status = data.get("conversion_loom_status")
                    conversion_rds_status = data.get("conversion_rds_status")
                    conversion_anndata_status = data.get("conversion_anndata_status")
                    self.assertIn(data.get("conversion_cxg_status"), expected_conversion_statuses)
                    if conversion_cxg_status == "FAILED":
                        self.fail(f"CXG CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_uuid}")
                    if conversion_loom_status == "FAILED":
                        self.fail(f"Loom CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_uuid}")
                    if conversion_rds_status == "FAILED":
                        self.fail(f"RDS CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_uuid}")
                    if conversion_anndata_status == "FAILED":
                        self.fail(f"Anndata CONVERSION FAILED. Status: {data}, Check logs for dataset: {dataset_uuid}")
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

        with self.subTest("test non owner cant retrieve status"):
            no_auth_headers = {"Content-Type": "application/json"}
            res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_uuid}/status", headers=no_auth_headers)
            with self.assertRaises(HTTPError):
                res.raise_for_status()

        with self.subTest("Test dataset deletion"):
            res = requests.delete(f"{self.api}/dp/v1/datasets/{dataset_uuid}", headers=headers)
            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.accepted)

            # Check that the dataset is gone
            res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_uuid}/status", headers=headers)
            self.assertEqual(res.status_code, requests.codes.forbidden)
