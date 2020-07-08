import unittest
from furl import furl
from datetime import datetime
from tests.unit.backend.chalice.api_server import BaseAPITest
import json


class TestProject(BaseAPITest, unittest.TestCase):
    def test__GET__project__OK(self):
        with self.subTest("No Parameters"):
            response = self.app.get("/v1/project", headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body["projects"][0]["id"], "test_project_id")

        with self.subTest("With from_date"):
            now = int(datetime.now().timestamp())
            test_url = furl(path="/v1/project", query_params={"from_date": now})
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual([], "test_project_id")
