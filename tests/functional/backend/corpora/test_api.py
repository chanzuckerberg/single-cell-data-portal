import json
import os
import unittest

import requests

from backend.corpora.common.utils.aws_secret import AwsSecret

if not os.getenv("DEPLOYMENT_STAGE"):  # noqa
    os.environ["DEPLOYMENT_STAGE"] = "dev"  # noqa

API_URL = {
    "dev": "https://api.dev.corpora.cziscience.com",
}


class TestApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.auth0_secret = json.loads(AwsSecret("corpora/test/auth0-secret").value)
        cls.auth0_secret["audience"] = API_URL.get(os.getenv("DEPLOYMENT_STAGE"))
        access_token = cls.get_auth_token()["access_token"]
        cls.auth_header = {"Authorization": f"bearer {access_token}"}

    def setUp(self):
        self.api = API_URL.get(os.environ["DEPLOYMENT_STAGE"])
        self.test_project_id = "005d611a-14d5-4fbf-846e-571a1f874f70"
        self.test_file_id = "7c93775542b056e048aa474535b8e5c2"
        self.bad_project_id = "DNE"
        self.bad_file_id = "DNE"

    def test_root_route(self):
        res = requests.get(f"{self.api}/")

        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)

    def test_get_projects(self):
        res = requests.get(f"{self.api}/projects")

        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)
        data = json.loads(res.content)

        for project in data:
            self.assertIs(type(project["id"]), str)
            self.assertIs(type(project["title"]), str)
            self.assertIs(type(project["assays"]), list)
            self.assertIs(type(project["organs"]), list)
            self.assertIs(type(project["species"]), list)
            self.assertIs(type(project["cell_count"]), int)

    @classmethod
    def get_auth_token(cls) -> dict:
        return requests.post(
            "https://czi-single-cell.auth0.com/oauth/token",
            json=cls.auth0_secret,
            headers={"content-type": "application/json"},
        ).json()
