import json
import os
import unittest

import requests
from dcplib.aws_secret import AwsSecret

if not os.getenv("DEPLOYMENT_STAGE"):  # noqa
    os.environ["DEPLOYMENT_STAGE"] = "test"  # noqa

API_URL = {
    "test": "https://browser-api-test.dev.single-cell.czi.technology",
    "dev": "https://browser-api.dev.single-cell.czi.technology",
}


class TestApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.auth0_secret = json.loads(
            AwsSecret(f"dcp/backend/browser/test/auth0-secret").value)
        cls.auth0_secret["audience"] = f"https://browser-api.{os.getenv('DEPLOYMENT_STAGE')}.single-cell.czi.technology"
        access_token = cls.get_auth_token()['access_token']
        cls.auth_header = {"Authorization": f"bearer {access_token}"}

    def setUp(self):
        self.api = API_URL[os.environ["DEPLOYMENT_STAGE"]]
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

    def test_get_project(self):
        with self.subTest("Project exists"):
            res = requests.get(f"{self.api}/projects/{self.test_project_id}")

            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.ok)

            test_project = json.loads(res.content)

            self.assertEqual(test_project["id"], "005d611a-14d5-4fbf-846e-571a1f874f70")
            self.assertEqual(test_project["title"], "HPSI human cerebral organoids")
            self.assertEqual(test_project["assays"], ["10X v2 sequencing"])
            self.assertEqual(test_project["organs"], ["skin of body", "brain", "stem cell"])
            self.assertEqual(test_project["species"], ["Homo sapiens"])
            self.assertEqual(test_project["cell_count"], 19900)

        with self.subTest("Project not found"):
            res = requests.get(f"{self.api}/projects/{self.bad_project_id}")
            self.assertEqual(res.status_code, requests.codes.not_found)

    def test_get_project_files(self):
        with self.subTest("Project exists"):
            res = requests.get(f"{self.api}/projects/{self.test_project_id}/files", headers=self.auth_header)

            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.ok)

            test_files = json.loads(res.content)
            self.assertEqual(len(test_files), 1)
            self.assertIs(type(test_files[0]["id"]), str)
            self.assertIs(type(test_files[0]["filename"]), str)
            self.assertIs(type(test_files[0]["file_format"]), str)
            self.assertIs(type(test_files[0]["file_type"]), str)
            self.assertIs(type(test_files[0]["file_size"]), int)

        with self.subTest("Project not found"):
            res = requests.get(f"{self.api}/projects/{self.bad_project_id}/files", headers=self.auth_header)
            self.assertEqual(res.status_code, requests.codes.not_found)

        with self.subTest("Not Authorized"):
            res = requests.get(f"{self.api}/projects/{self.test_project_id}/files")
            self.assertEqual(res.status_code, requests.codes.unauthorized)

    def test_get_file(self):
        with self.subTest("File exists"):
            res = requests.get(f"{self.api}/files/{self.test_file_id}", headers=self.auth_header)

            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.ok)

        with self.subTest("File not found"):
            res = requests.get(f"{self.api}/files/{self.bad_file_id}", headers=self.auth_header)
            self.assertEqual(res.status_code, requests.codes.not_found)

        with self.subTest("Not Authorized"):
            res = requests.get(f"{self.api}/files/{self.test_file_id}",headers=self.auth_header)
            self.assertEqual(res.status_code, requests.codes.unauthorized)

    @classmethod
    def get_auth_token(cls) -> dict:
        return requests.post(
            "https://czi-single-cell.auth0.com/oauth/token",
            json=cls.auth0_secret,
            headers={"content-type": "application/json"},
        ).json()
