import base64
import json
import os
import unittest
import requests
from requests import HTTPError

from backend.corpora.common.corpora_config import CorporaAuthConfig

API_URL = {"dev": "https://api.dev.corpora.cziscience.com", "test": "http://localhost:5000"}


class TestApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.deployment_stage = os.environ["DEPLOYMENT_STAGE"]

    def setUp(self):
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
                username="user@example.com",
                password=config.test_account_password,
                audience="cellxgene.dev.single-cell.czi.technology/",
                scope="openid profile email offline",
                client_id=config.client_id,
                client_secret=config.client_secret,
            ),
        )

        id_token = response.json()["id_token"]
        token = {"id_token": id_token}
        cls.cookie = base64.b64encode(json.dumps(dict(token)).encode("utf-8")).decode()

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
        self.assertEqual(res.status_code, requests.codes.created)
        data = json.loads(res.content)
        self.assertIn("collection_uuid", data)

        collection_uuid = data["collection_uuid"]
        # check created collection returns as private
        res = requests.get(f"{self.api}/dp/v1/collections", headers=headers)
        data = json.loads(res.content)
        private_collection_uuids = []
        for collection in data["collections"]:
            if collection["visibility"] == "PRIVATE":
                private_collection_uuids.append(collection["id"])

        self.assertIn(collection_uuid, private_collection_uuids)

        # make collection public
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

        # check collection available to everyone
        no_auth_headers = {"Content-Type": "application/json"}
        res = requests.get(f"{self.api}/dp/v1/collections", headers=no_auth_headers)
        data = json.loads(res.content)
        collection_uuids = [x["id"] for x in data["collections"]]
        self.assertIn(collection_uuid, collection_uuids)

        # TODO @madison -- Add once we can delete collections
        # delete collection
        # res = requests.delete(f"{self.api}/dp/v1/collections/{collection_uuid}", headers=headers)
        # res.raise_for_status()
        # self.assertEqual(res.status_code, requests.codes.accepted)
        #
        # # check collection gone
        # res = requests.get(f"{self.api}/dp/v1/collections", headers=no_auth_headers)
        # data = json.loads(res.content)
        # collection_uuids = [x['id'] for x in data['collections']]
        # self.assertNotIn(collection_uuid, collection_uuids)

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
        self.assertEqual(res.status_code, requests.codes.created)
        data = json.loads(res.content)
        self.assertIn("collection_uuid", data)

        collection_uuid = data["collection_uuid"]

        body = {"url": "https://www.dropbox.com/s/ib9pth7jr5fvaa8/7MB.h5ad?dl=0"}

        res = requests.post(
            f"{self.api}/dp/v1/collections/{collection_uuid}/upload-links", data=json.dumps(body), headers=headers
        )
        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.accepted)

        dataset_uuid = json.loads(res.content)["dataset_uuid"]
        res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_uuid}/status", headers=headers)
        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)
        data = json.loads(res.content)
        # TODO @madison update if we speed up backend
        self.assertEqual(data["upload_status"], "WAITING")

        # Check non_auth user cant check status
        no_auth_headers = {"Content-Type": "application/json"}
        res = requests.get(f"{self.api}/dp/v1/datasets/{dataset_uuid}/status", headers=no_auth_headers)
        with self.assertRaises(HTTPError):
            res.raise_for_status()
