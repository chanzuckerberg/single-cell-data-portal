from tests.unit.backend.chalice.api_server.base_api_test import BaseAPITest
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestAPI(BaseAPITest, DataPortalTestCase):
    def test_smoke(self):
        """ If this fails then the server does not work """
        response = self.app.get("/")
        response.raise_for_status()
        self.assertIn("X-AWS-REQUEST-ID", response.headers.keys())
