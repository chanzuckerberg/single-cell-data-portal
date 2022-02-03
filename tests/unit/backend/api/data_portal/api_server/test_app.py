from unit.backend.api.data_portal.api_server.base_api_test import BaseAPITest


class TestAPI(BaseAPITest):
    def test_smoke(self):
        """ If this fails then the server does not work """
        response = self.app.get("/")
        self.assertEqual(200, response.status_code)
