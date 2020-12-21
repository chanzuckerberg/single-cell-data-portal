import json
import os
import unittest

import sys
from furl import furl

from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.corpora import CorporaTestCaseUsingMockAWS
from tests.unit.backend.chalice.api_server.mock_auth import MockOauthServer, get_auth_token


class TestDataset(BaseAPITest, CorporaTestCaseUsingMockAWS):
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

    def test__post_dataset_asset__OK(self):
        bucket = self.CORPORA_TEST_CONFIG["bucket_name"]
        s3_file_name = "test_s3_uri.h5ad"
        content = "Hello world!"
        self.create_s3_object(s3_file_name, bucket, content=content)

        expected_body = dict(dataset_id="test_dataset_id", file_name="test_filename", file_size=len(content))
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        response.raise_for_status()
        actual_body = json.loads(response.body)
        presign_url = actual_body.pop("presigned_url")
        self.assertIsNotNone(presign_url)
        self.assertEqual(expected_body, actual_body)

    @unittest.skipIf(os.getenv("IS_DOCKER_DEV"), "This is currently TODO in docker")
    def test__post_dataset_asset__file_SERVER_ERROR(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(500, response.status_code)
        body = json.loads(response.body)
        self.assertEqual("An internal server error has occurred. Please try again later.", body["detail"])

    def test__post_dataset_asset__dataset_NOT_FOUND(self):
        test_url = furl(path="/dp/v1/datasets/test_user_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.body)
        self.assertEqual("'dataset/test_user_id' not found.", body["detail"])
        print(body)

    def test__post_dataset_asset__asset_NOT_FOUND(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/fake_asset")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.body)
        self.assertEqual("'dataset/test_dataset_id/asset/fake_asset' not found.", body["detail"])

    def test__get_status__ok(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        response.raise_for_status()
        actual_body = json.loads(response.body)
        expected_body = {
            "conversion_anndata_status": "NA",
            "conversion_cxg_status": "NA",
            "conversion_loom_status": "NA",
            "conversion_rds_status": "NA",
            "dataset_id": "test_dataset_id",
            "id": "test_dataset_processing_status_id",
            "upload_progress": 0.4444444444444444,
            "upload_status": "UPLOADING",
            "validation_status": "NA",
        }
        self.assertEqual(expected_body, actual_body)

    def test__get_status__403(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id_not_owner/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(403, response.status_code)
