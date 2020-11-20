import json
import sys
import os
import unittest

from furl import furl

from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.chalice.api_server.mock_auth import MockOauthServer, get_auth_token


class TestCollectionUploadLink(BaseAPITest, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        BaseAPITest.setUpClass()
        cls.mock_oauth_server = MockOauthServer()
        cls.mock_oauth_server.start()
        assert cls.mock_oauth_server.server_okay

        cls.old_path = sys.path.copy()
        sys.path.insert(0, os.path.join(cls.corpora_api_dir, "chalicelib"))  # noqa
        from corpora.common.corpora_config import CorporaAuthConfig

        # Use the CorporaAuthConfig used by the chalice app
        cls.auth_config = CorporaAuthConfig()
        cls.auth_config._config["api_base_url"] = f"http://localhost:{cls.mock_oauth_server.port}"
        cls.auth_config._config["callback_base_url"] = "http://localhost:5000"
        cls.auth_config.update_defaults()

    @classmethod
    def tearDownClass(cls):
        cls.mock_oauth_server.terminate()
        sys.path = cls.old_path

    def test__link__Accepted(self):
        path = "/dp/v1/collections/test_collection_id/upload/link"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {
            'url': 'https://test_url.com'
        }

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        response.raise_for_status()
        actual_body = json.loads(response.body)
