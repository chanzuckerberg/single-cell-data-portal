from tests.unit.backend.corpora.api_server.base_api_test import BaseAPITest


class TestAPI(BaseAPITest):
    def test_smoke(self):
        """ If this fails then the server does not work """
        response = self.app.get("/")
        response.raise_for_status()
