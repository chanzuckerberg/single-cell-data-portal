import logging

from tests.unit.backend.layers.common.base_api_test import NewBaseTest


class TestAPI(NewBaseTest):
    def check_request_id(self, response):
        request_id = [header[1] for header in response.headers if header[0] == "X-Request-Id"][0]
        self.assertTrue(len(request_id) > 1)

    def test__api_index__returns_200(self):
        """If this fails then the server does not work"""
        response = self.app.get("/")
        self.assertEqual(200, response.status_code)
        self.check_request_id(response)

    def test__dp_api_swagger_ui__returns_200(self):
        response = self.app.get("/dp/ui/")
        self.assertEqual(200, response.status_code)

    def test__curation_api_swagger_ui__returns_200(self):
        response = self.app.get("/curation/ui/")
        self.assertEqual(200, response.status_code)

    def test__wmg_api_swagger_ui__returns_200(self):
        response = self.app.get("/wmg/ui/")
        self.assertEqual(200, response.status_code)

    def test__catchall_exception_handler(self):
        with self.assertLogs("corpora-api-test", level=logging.ERROR) as logs:
            response = self.app.get("/exception")
        self.assertTrue(len(logs.records) == 1)
        self.assertEqual(response.status_code, 500)
        self.check_request_id(response)
