import json
import unittest

from furl import furl

from backend.corpora.common.corpora_orm import ProjectStatus
from backend.corpora.common.entities import Project
from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.utils import BogusProjectParams
from datetime import datetime


class TestProject(BaseAPITest, unittest.TestCase):
    def test__list_project__ok(self):
        with self.subTest("No Parameters"):
            test_url = furl(path="/v1/project")
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertIn("test_project_id", [p["id"] for p in body["projects"]])

        creation_time = 0

        test_project = Project.create(
            **BogusProjectParams.get(status=ProjectStatus.LIVE.name), created_at=datetime.fromtimestamp(creation_time)
        )
        test_id = test_project.id
        future_time = int(datetime.fromtimestamp(10).timestamp())
        with self.subTest("With to_date"):
            test_url = furl(path="/v1/project", query_params={"to_date": future_time})
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["projects"][0]["id"], test_id)
            self.assertEqual(body["projects"][0]["created_at"], creation_time)
            self.assertEqual(body["to_date"], future_time)

        with self.subTest("With from_date"):
            test_url = furl(path="/v1/project", query_params={"from_date": future_time})
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertIn("test_project_id", [p["id"] for p in body["projects"]])
            self.assertEqual(body["from_date"], future_time)

    def test__get_project_uuid__ok(self):
        with self.subTest("Exists"):
            test_url = furl(path="/v1/project/test_project_id")
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["id"], "test_project_id")
