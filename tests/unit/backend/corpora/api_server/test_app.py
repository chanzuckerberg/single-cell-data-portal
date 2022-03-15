from tests.unit.backend.corpora.api_server.base_api_test import BaseAPITest


class TestAPI(BaseAPITest):
    def test_smoke_dataportal(self):
        """If this fails then the server does not work"""
        response = self.app.get("/ui/")
        self.assertEqual(200, response.status_code)

    def test_smoke_curator(self):
        response = self.app.get("/curator/v1/ui/")
        self.assertEqual(200, response.status_code)
