import json
import os
import unittest

import requests

API_URL = {
    "dev": "https://api.dev.corpora.cziscience.com",
    "test": "http://localhost:3000"
}


class TestApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.deployment_stage = os.getenv("DEPLOYMENT_STAGE", "test")

    def setUp(self):
        self.api = API_URL.get(self.deployment_stage)
        self.test_project_id = "005d611a-14d5-4fbf-846e-571a1f874f70"
        self.test_file_id = "7c93775542b056e048aa474535b8e5c2"
        self.bad_project_id = "DNE"
        self.bad_file_id = "DNE"

    def test_root_route(self):
        res = requests.get(f"{self.api}/")

        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)

    def test_get_projects(self):
        res = requests.get(f"{self.api}/v1/project")

        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)
        data = json.loads(res.content)

        for project in data['projects']:
            self.assertIsInstance(project["id"], str)
            self.assertIsInstance(project["created_at"], float)
