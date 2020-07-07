import unittest

from tests.unit.backend.chalice.api_server import BaseAPITest


class TestProject(BaseAPITest, unittest.TestCase):
    def test__GET__project__OK(self):
        response = self.app.get("/v1/project", headers=dict(host='localhost'))
        response.raise_for_status()
