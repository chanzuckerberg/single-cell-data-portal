import functools
import json
import os
import unittest

import requests
from chalice.cli import CLIFactory
from chalice.local import LocalGateway, LocalGatewayException


class TestApi(unittest.TestCase):
    def setUp(self):
        self.app = ChaliceTestHarness()
        self.test_project_id = "005d611a-14d5-4fbf-846e-571a1f874f70"
        self.test_file_id = "7c93775542b056e048aa474535b8e5c2"
        self.bad_project_id = "DNE"
        self.bad_file_id = "DNE"

    def test_root_route(self):
        res = self.app.get("/")

        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)

    def test_get_projects(self):
        res = self.app.get("/projects")

        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)

        data = json.loads(res.body)

        for project in data:
            self.assertIs(type(project["id"]), str)
            self.assertIs(type(project["title"]), str)
            self.assertIs(type(project["assays"]), list)
            self.assertIs(type(project["organs"]), list)
            self.assertIs(type(project["species"]), list)
            self.assertIs(type(project["cell_count"]), int)

    def test_get_project(self):
        with self.subTest("Project exists"):
            res = self.app.get(f"/projects/{self.test_project_id}")

            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.ok)

            test_project = json.loads(res.body)

            self.assertEqual(test_project["id"], "005d611a-14d5-4fbf-846e-571a1f874f70")
            self.assertEqual(test_project["title"], "HPSI human cerebral organoids")
            self.assertEqual(test_project["assays"], ["10X v2 sequencing"])
            self.assertEqual(test_project["organs"], ["skin of body", "brain", "stem cell"])
            self.assertEqual(test_project["species"], ["Homo sapiens"])
            self.assertEqual(test_project["cell_count"], 19900)

        with self.subTest("Project not found"):
            res = self.app.get(f"/projects/{self.bad_project_id}")
            self.assertEqual(res.status_code, requests.codes.not_found)

    def test_get_project_files(self):
        with self.subTest("Project exists"):
            res = self.app.get(f"/projects/{self.test_project_id}/files")

            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.ok)

            test_files = json.loads(res.body)
            self.assertEqual(len(test_files), 1)
            self.assertIs(type(test_files[0]["id"]), str)
            self.assertIs(type(test_files[0]["filename"]), str)
            self.assertIs(type(test_files[0]["file_format"]), str)
            self.assertIs(type(test_files[0]["file_type"]), str)
            self.assertIs(type(test_files[0]["file_size"]), int)

        with self.subTest("Project not found"):
            res = self.app.get(f"/projects/{self.bad_project_id}/files")
            self.assertEqual(res.status_code, requests.codes.not_found)

    def test_get_file(self):
        with self.subTest("File exists"):
            res = self.app.get(f"/files/{self.test_file_id}")

            res.raise_for_status()
            self.assertEqual(res.status_code, requests.codes.ok)

        with self.subTest("File not found"):
            res = self.app.get(f"/files/{self.bad_file_id}")
            self.assertEqual(res.status_code, requests.codes.not_found)


class ChaliceTestHarness:
    def __init__(self):
        project_dir = os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "..", "..", "dcp_prototype", "backend", "browser", "api"
        )
        config = CLIFactory(project_dir=project_dir).create_config_obj(chalice_stage_name="dev")
        self._chalice_app = config.chalice_app
        self._gateway = LocalGateway(self._chalice_app, config)

    @functools.lru_cache(maxsize=64, typed=False)
    def __getattr__(self, method):
        return functools.partial(self.request, method=method.upper())

    def request(self, path, headers={}, body={}, method="GET"):
        resp_obj = requests.Response()
        try:
            response = self._gateway.handle_request(method, path, headers, body)
        except LocalGatewayException as error:
            resp_obj.status_code = error.CODE
            resp_obj.headers = error.headers
            resp_obj.body = error.body
        else:
            resp_obj.status_code = response["statusCode"]
            resp_obj.headers = response["headers"]
            resp_obj.body = response["body"]
        resp_obj.headers["Content-Length"] = str(len(body))
        return resp_obj


if __name__ == "__main__":
    unittest.main()
