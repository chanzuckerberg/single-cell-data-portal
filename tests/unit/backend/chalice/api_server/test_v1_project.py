import unittest
import furl
from tests.unit.backend.chalice.api_server import BaseAPITest
import json


class TestProject(BaseAPITest, unittest.TestCase):

    def test__GET__project__OK(self):
        with self.subTest("No Parameters"):
            response = self.app.get("/v1/project", headers=dict(host="localhost"))
            response.raise_for_status()
            body = json.loads(response.body)
            self.assertEqual(body['projects'][0]['id'], 'test_project_id')
