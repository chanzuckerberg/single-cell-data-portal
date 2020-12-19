import json
import unittest
from mock import patch
from furl import furl

from backend.corpora.common.utils.math_utils import GB
from tests.unit.backend.chalice.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.chalice.api_server.mock_auth import get_auth_token
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup, fixture_file_path


class TestCollectionUploadLink(BaseAuthAPITest, unittest.TestCase):
    def setUp(self):
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        self.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"

    @patch("corpora.common.upload_sfn.start_upload_sfn")
    def test__link__202(self, mocked):
        with EnvironmentSetup({"CORPORA_CONFIG": fixture_file_path("bogo_config.js")}):
            path = "/dp/v1/collections/test_collection_id/upload-links"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            body = {"url": self.good_link}

            test_url = furl(path=path)
            response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.assertIn("dataset_uuid", actual_body.keys())

    def test__link_no_auth__401(self):
        path = "/dp/v1/collections/test_collection_id/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(401, response.status_code)

    @patch("corpora.common.utils.dropbox.get_file_info", return_value={"size": 1, "name": "file.h5ad"})
    def test__link_not_owner__403(self, mock_get_file_info):
        path = "/dp/v1/collections/test_collection_id_not_owner/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    def test__bad_link__400(self):
        path = "/dp/v1/collections/test_collection_id/upload-links"

        with self.subTest("Unsupported Provider"):
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            body = {"url": "https://test_url.com"}
            test_url = furl(path=path)
            response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
            self.assertEqual(400, response.status_code)
            self.assertEqual("The dropbox shared link is invalid.", json.loads(response.body)["detail"])

        with self.subTest("Bad Dropbox link"):
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            body = {"url": self.dummy_link}
            test_url = furl(path=path)
            response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
            self.assertEqual(400, response.status_code)
            self.assertEqual("The URL provided causes an error with Dropbox.", json.loads(response.body)["detail"])

    @patch("corpora.common.utils.dropbox.get_file_info", return_value={"size": 1, "name": "file.txt"})
    def test__unsupported_format__400(self, mock_func):
        path = "/dp/v1/collections/test_collection_id/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(400, response.status_code)

    @patch("corpora.common.utils.dropbox.get_file_info", return_value={"size": 31 * GB, "name": "file.txt"})
    def test__oversized__413(self, mock_func):
        path = "/dp/v1/collections/test_collection_id/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(413, response.status_code)

    def test__link_fake_collection__403(self):
        path = "/dp/v1/collections/fake/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.good_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    def test_link_live_collection__403(self):
        path = "/dp/v1/collections/test_collection_id_public/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.good_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)
