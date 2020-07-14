import json
import unittest

from furl import furl

from backend.corpora.common.entities import Project
from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.utils import ProjectParams


class TestProject(BaseAPITest, unittest.TestCase):
    def test__get_project__ok(self):
        test_id = Project.create(**ProjectParams.get(), created_at=10, updated_at=10).id
        with self.subTest("No Parameters"):
            test_url = furl(path="/v1/project")
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["projects"][0]["id"], test_id)

        now = 15
        with self.subTest("With to_date"):
            test_url = furl(path="/v1/project", query_params={"to_date": now})
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["projects"][0]["id"], test_id)

        with self.subTest("With from_date"):
            test_url = furl(path="/v1/project", query_params={"from_date": now})
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["projects"], [])

    def test__get_project_uuid__ok(self):
        with self.subTest("Exists"):
            test_url = furl(path="/v1/project/test_project_id")
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["id"], "test_project_id")
