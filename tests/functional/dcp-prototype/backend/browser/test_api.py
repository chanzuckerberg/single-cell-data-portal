import json
import os
import unittest

import requests

API_URL = {
    "test": "https://c5eyvsi657.execute-api.us-east-1.amazonaws.com/test",
    "dev": "https://u5j1wa9u5i.execute-api.us-east-1.amazonaws.com/dev"
}


class TestApi(unittest.TestCase):
    def setUp(self):
        self.api = API_URL[os.environ['DEPLOYMENT_STAGE']]
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
            res = requests.get(f"{self.api}/projects/{self.test_project_id}/files")

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
            res = requests.get(f"{self.api}/projects/{self.bad_project_id}/files")
            self.assertEqual(res.status_code, requests.codes.not_found)

    def test_get_file(self):
        with self.subTest("File exists"):
            res = requests.get(f"{self.api}/files/{self.test_file_id}")

            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.ok)

        with self.subTest("File not found"):
            res = requests.get(f"{self.api}/files/{self.bad_file_id}")
            self.assertEqual(res.status_code, requests.codes.not_found)
