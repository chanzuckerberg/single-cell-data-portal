import json
import os
import unittest

import requests

API_URL = {"dev": "https://api.dev.corpora.cziscience.com", "test": "http://localhost:5000"}


class TestApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.deployment_stage = os.environ["DEPLOYMENT_STAGE"]

    def setUp(self):
        self.api = API_URL.get(self.deployment_stage)
        self.test_collection_id = "005d611a-14d5-4fbf-846e-571a1f874f70"
        self.test_file_id = "7c93775542b056e048aa474535b8e5c2"
        self.bad_collection_id = "DNE"
        self.bad_file_id = "DNE"

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

    def test_create_and_retrieve_collection(self):
        data = json.dumps(dict(
            name="test collection",
            description="This is a test collection",
            contact_name="person human",
            contact_email="person@human.com",
            data_submission_policy_version="0.0.1",
            links=[
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"}]
        ))
        res = requests.post(
            f"{self.api}/dp/v1/collections",
            headers={"host": "localhost", 'Content-Type': "application/json"},
            data=data
        )
        import pdb
        pdb.set_trace()
        self.assertEqual(res.status_code, requests.codes.created)
        import pdb
        pdb.set_trace()
