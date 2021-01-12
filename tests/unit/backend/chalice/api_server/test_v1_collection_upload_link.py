import json
import sys
import os
from mock import patch
from furl import furl

from backend.corpora.common.corpora_orm import UploadStatus
from backend.corpora.common.utils.math_utils import GB
from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.chalice.api_server.generate_data_mixin import GenerateDataMixin
from tests.unit.backend.chalice.api_server.mock_auth import MockOauthServer, get_auth_token
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup, fixture_file_path
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestCollectionUploadLink(BaseAPITest, DataPortalTestCase, GenerateDataMixin):
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
        super().setUp()
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

    def test__cancel_dataset_download__ok(self):
        # Test pre upload
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)
        self.assertEqual(json.loads(response.body)["upload_status"], "CANCEL_PENDING")

        # Test while uploading
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        dataset = self.generate_dataset(processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)
        self.assertEqual(json.loads(response.body)["upload_status"], "CANCEL_PENDING")

    def test__cancel_dataset_download__dataset_does_not_exist(self):
        test_url = "/dp/v1/datasets/missing_dataset_id"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__cancel_dataset_download__already_uploaded(self):
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 0.0}
        dataset = self.generate_dataset(processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 405)
        self.assertEqual(
            json.loads(response.body)["detail"], f"'dataset/{dataset.id}' upload is complete and can not be cancelled."
        )

    def test__cancel_dataset_download__user_not_collection_owner(self):
        collection = self.generate_collection(owner="someone_else")
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__cancel_dataset_download__user_not_logged_in(self):
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 401)
