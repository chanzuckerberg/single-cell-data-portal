from tests.unit.backend.corpora.api_server.base_api_test import BaseAPITest


class TestAPI(BaseAPITest):
    def test__api_index__returns_200(self):
        """If this fails then the server does not work"""
        response = self.app.get("/")
        self.assertEqual(200, response.status_code)
        request_id = response.headers[2]
        self.assertTrue(request_id[0] == "X-Request-Id")
        self.assertTrue(len(request_id[1]) > 0)

    def test__dp_api_swagger_ui__returns_200(self):
        response = self.app.get("/dp/ui/")
        self.assertEqual(200, response.status_code)

    def test__curation_api_swagger_ui__returns_200(self):
        response = self.app.get("/curation/ui/")
        self.assertEqual(200, response.status_code)

    def test__wmg_api_swagger_ui__returns_200(self):
        response = self.app.get("/wmg/ui/")
        self.assertEqual(200, response.status_code)
