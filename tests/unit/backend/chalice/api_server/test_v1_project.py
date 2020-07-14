import json
import unittest

from furl import furl

from backend.corpora.common.entities import Project
from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.utils import ProjectParams
from datetime import datetime

class TestProject(BaseAPITest, unittest.TestCase):
    def test__list_project__ok(self):
        with self.subTest("No Parameters"):
            test_url = furl(path="/v1/project")
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertIn("test_project_id", [p["id"] for p in body["projects"]])

        creation_time = 10
        test_id = Project.create(**ProjectParams.get(), created_at=datetime.utcfromtimestamp(creation_time)).id

        now = 15
        with self.subTest("With to_date"):
            test_url = furl(path="/v1/project", query_params={"to_date": now})
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["projects"][0]["id"], test_id)
            self.assertEqual(body["projects"][0]["created_at"], creation_time)
            self.assertEqual(body["to_date"], now)


        with self.subTest("With from_date"):
            test_url = furl(path="/v1/project", query_params={"from_date": now})
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["projects"], [])
            self.assertEqual(body["from_date"], now)

    def test__get_project_uuid__ok(self):
        with self.subTest("Exists"):
            test_url = furl(path="/v1/project/test_project_id")
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["id"], "test_project_id")
