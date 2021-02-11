import json
import os
import unittest

from furl import furl

from backend.corpora.common.corpora_orm import UploadStatus, CollectionVisibility
from backend.corpora.common.utils.db_session import processing_status_updater
from tests.unit.backend.chalice.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.chalice.api_server.mock_auth import get_auth_token
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestDataset(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        CorporaTestCaseUsingMockAWS.setUp(self)

    def tearDown(self):
        CorporaTestCaseUsingMockAWS.tearDown(self)

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
            "processing_status": "PENDING",
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

    def test__minimal_status__ok(self):
        dataset = self.generate_dataset(
            self.session, processing_status={"upload_status": "WAITING", "upload_progress": 0.0}
        )
        processing_status_id = dataset.processing_status.id
        test_url = furl(path=f"/dp/v1/datasets/{dataset.id}/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        response.raise_for_status()
        actual_body = json.loads(response.body)
        expected_body = {
            "dataset_id": dataset.id,
            "id": processing_status_id,
            "upload_progress": 0.0,
            "upload_status": "WAITING",
        }
        self.assertEqual(expected_body, actual_body)

        for status in UploadStatus:
            processing_status = {"upload_status": status, "upload_progress": 0.0}
            processing_status_updater(self.session, processing_status_id, processing_status)
            response = self.app.get(test_url.url, headers=headers)
            self.assertEqual(json.loads(response.body)["upload_status"], status.name)

    def test__cancel_dataset_download__ok(self):
        # Test pre upload
        collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # Test while uploading
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        dataset = self.generate_dataset(collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

    def test__cancel_dataset_download__dataset_does_not_exist(self):
        test_url = "/dp/v1/datasets/missing_dataset_id"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__delete_uploaded_dataset__ok(self):
        collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 0.0}
        dataset = self.generate_dataset(collection=collection, processing_status=processing_status)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        # check dataset in collection
        collection_url = furl(path=f"/dp/v1/collections/{collection.id}")
        collection_url.add(query_params=dict(visibility=CollectionVisibility.PRIVATE.name))
        response = self.app.get(collection_url.url, headers=headers)
        response.raise_for_status()
        body = json.loads(response.body)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset.id, dataset_ids)

        # delete dataset
        test_url = f"/dp/v1/datasets/{dataset.id}"
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

        # check dataset no longer returned in collection
        response = self.app.get(collection_url.url, headers=headers)
        response.raise_for_status()
        body = json.loads(response.body)

        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertNotIn(dataset.id, dataset_ids)

    def test__delete_dataset_can_safetly_be_called_twice(self):
        collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        dataset = self.generate_dataset(collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # delete again
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

    def test__get_deleted_dataset_status__returns_403(self):
        collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id")
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 0.0}
        dataset = self.generate_dataset(collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}/status"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url, headers=headers)
        response.raise_for_status()
        actual_body = json.loads(response.body)
        expected_body = {
            "dataset_id": dataset.id,
            "id": dataset.processing_status.id,
            "upload_progress": 0.0,
            "upload_status": "UPLOADED",
        }
        self.assertEqual(expected_body, actual_body)

        # delete the dataset
        self.app.delete(f"/dp/v1/datasets/{dataset.id}", headers=headers)

        response = self.app.get(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__delete_public_dataset_returns__405(self):
        collection = self.generate_collection(visibility=CollectionVisibility.PUBLIC.name)
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 0.0}
        dataset = self.generate_dataset(collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 405)
        self.assertEqual(json.loads(response.body), "Can not delete a public dataset")

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
