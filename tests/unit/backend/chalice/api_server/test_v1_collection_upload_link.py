import json
import sys
import os
import unittest
from mock import patch
from furl import furl

from backend.corpora.common.utils.math_utils import GB
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

    def setUp(self):
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        self.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"

    @patch("corpora.common.upload_sfn.start_upload_sfn")
    def test__link__accepted(self, mocked):
        path = "/dp/v1/collections/test_collection_id/upload/link"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.good_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        response.raise_for_status()
        actual_body = json.loads(response.body)
        self.assertIn("dataset_uuid", actual_body.keys())

    def test__link_no_auth__401(self):
        path = "/dp/v1/collections/test_collection_id/upload/link"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(401, response.status_code)

    @patch("corpora.common.utils.dropbox.get_file_info", return_value={"size": 1, "name": "file.h5ad"})
    def test__link_not_owner__403(self, mock_get_file_info):
        path = "/dp/v1/collections/test_collection_id_not_owner/upload/link"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    def test__bad_link__400(self):
        path = "/dp/v1/collections/test_collection_id/upload/link"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": "https://test_url.com"}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(400, response.status_code)

    @patch("corpora.common.utils.dropbox.get_file_info", return_value={"size": 1, "name": "file.txt"})
    def test__unsupported_format__400(self, mock_func):
        path = "/dp/v1/collections/test_collection_id/upload/link"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(400, response.status_code)

    @patch("corpora.common.utils.dropbox.get_file_info", return_value={"size": 31 * GB, "name": "file.txt"})
    def test__oversized__413(self, mock_func):
        path = "/dp/v1/collections/test_collection_id/upload/link"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(413, response.status_code)
