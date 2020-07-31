import unittest

from tests.unit.backend.chalice.api_server import BaseAPITest


class TestAPI(BaseAPITest, unittest.TestCase):
    def test_smoke(self):
        """ If this fails then the server does not work """
        response = self.app.get("/")
        response.raise_for_status()
        self.assertIn("X-AWS-REQUEST-ID", response.headers.keys())
