import json
import unittest

import boto3
from furl import furl
from moto import mock_s3

from tests.unit.backend.chalice.api_server import BaseAPITest


class TestDataset(BaseAPITest, unittest.TestCase):
    @mock_s3
    def test__post_dataset_artifact__OK(self):
        bucket = "corpora-data-test"
        region = "us-west-2"
        s3_file_name = "test_s3_uri.h5ad"
        file_body = "Hello world!"
        conn = boto3.resource("s3", region_name=region)
        conn.create_bucket(Bucket=bucket)
        s3 = boto3.client("s3", region_name=region)
        s3.put_object(Bucket=bucket, Key=s3_file_name, Body=file_body)

        expected_body = dict(
            dataset_id="test_dataset_id", file_name="test_filename", file_type="H5AD", file_size=len(file_body)
        )
        test_url = furl(path="/v1/dataset/test_dataset_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        response.raise_for_status()
        actual_body = json.loads(response.body)
        presign_url = actual_body.pop("presigned_url")
        self.assertIsNotNone(presign_url)
        self.assertEqual(expected_body, actual_body)

    @mock_s3
    def test__post_dataset_artifact__file_NOT_FOUND(self):
        bucket = "corpora-data-test"
        region = "us-west-2"
        conn = boto3.resource("s3", region_name=region)
        conn.create_bucket(Bucket=bucket)

        test_url = furl(path="/v1/dataset/test_dataset_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.body)
        self.assertEqual("'dataset/test_dataset_id/asset/test_dataset_artifact_id' not found.", body['detail'])

    def test__post_dataset_artifact__dataset_NOT_FOUND(self):
        test_url = furl(path="/v1/dataset/fake_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.body)
        self.assertEqual("'dataset/fake_id' not found.", body['detail'])
        print(body)

    def test__post_dataset_artifact__asset_NOT_FOUND(self):
        test_url = furl(path="/v1/dataset/test_dataset_id/asset/fake_asset")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.body)
        self.assertEqual("'dataset/test_dataset_id/asset/fake_asset' not found.", body['detail'])
