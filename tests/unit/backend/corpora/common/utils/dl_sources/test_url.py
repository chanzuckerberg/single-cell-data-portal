import os
import unittest
from unittest import TestCase
from urllib.parse import urlparse

import boto3
import requests
from moto import mock_s3

from backend.common.utils.dl_sources.url import from_url, RegisteredSources, S3URL, DropBoxURL, S3URI
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestRegisteredSources(unittest.TestCase):
    def test_register_neg(self):
        with self.assertRaises(TypeError):
            RegisteredSources.add(123)

    def test_register_pos(self):
        self.assertIn(S3URL, RegisteredSources.get())
        RegisteredSources.remove(S3URL)
        self.assertNotIn(S3URL, RegisteredSources.get())
        RegisteredSources.add(S3URL)
        self.assertIn(S3URL, RegisteredSources.get())

    def test_from_url(self):
        RegisteredSources.add(S3URL)
        self.assertIsInstance(from_url("https://somebucket.s3.amazonaws.com"), S3URL)

        RegisteredSources.add(DropBoxURL)
        self.assertIsInstance(from_url("https://www.dropbox.com/s/abcd1234/abcd1234?dl=0"), DropBoxURL)
        self.assertIsNone(from_url("http://example.com"))

        RegisteredSources.add(S3URI)
        self.assertIsInstance(from_url("s3://bucket/key"), S3URI)


class Test_S3URL(CorporaTestCaseUsingMockAWS):
    @staticmethod
    def get_presigned_url(dataset_id, asset_id):
        response = requests.post(
            f"https://api.cellxgene.staging.single-cell.czi.technology/dp/v1/datasets/" f"{dataset_id}/asset/{asset_id}"
        )
        response.raise_for_status()
        return response.json()["presigned_url"]

    def test_S3URL(self):
        # Reset the S3URL class back to normal after test
        self.addCleanup(setattr, S3URL, "_netloc", S3URL._netloc)
        self.addCleanup(setattr, S3URL, "_scheme", S3URL._scheme)

        # Setup S3URL to work with localstack
        parsed_url = urlparse(os.getenv("BOTO_ENDPOINT_URL"))
        S3URL._netloc = parsed_url.netloc
        S3URL._scheme = parsed_url.scheme

        # Store a file in the s3 bucket
        bucket = self.CORPORA_TEST_CONFIG["bucket_name"]
        s3_file_name = "test_s3_uri.h5ad"
        content = "Hello world!"
        self.create_s3_object(s3_file_name, bucket, content=content)

        # Generate presigned URL for the s3 object.
        s3 = boto3.client("s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"))
        presigned_url = s3.generate_presigned_url(
            "get_object", Params={"Bucket": self.corpora_config.bucket_name, "Key": s3_file_name}, ExpiresIn=3000
        )

        # Test S3URL
        url = S3URL.validate(presigned_url)
        self.assertIsInstance(url, S3URL)
        file_info = url.file_info()
        self.assertEqual(file_info["size"], 12)
        self.assertTrue(file_info["name"].endswith(s3_file_name))


class TestS3URI(TestCase):
    def test__validate_with_valid_s3_uri__ok(self):
        s3_uri = S3URI.validate("s3://bucket/key")

        self.assertEqual("s3", s3_uri.scheme)
        self.assertEqual("bucket", s3_uri.netloc)
        self.assertEqual("/key", s3_uri.path)

    def test__validate_with_invalid_s3_uri__returns_none(self):
        s3_uri = S3URI.validate("s4://bucket/key")
        self.assertIsNone(s3_uri)

        s3_uri = S3URI.validate("s3://bucket")
        self.assertIsNone(s3_uri)

    @mock_s3
    def test__file_info__returns_file_size(self):
        s3 = boto3.client("s3")
        bucket_name = "bucket"
        key = "/key"
        content = "stuff"
        s3.create_bucket(
            Bucket=bucket_name, CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]}
        )
        s3.put_object(Bucket=bucket_name, Key=key, Body=content)

        s3_uri = from_url(f"s3://{bucket_name}{key}")
        info = s3_uri.file_info()

        self.assertEqual(key, info["name"])
        self.assertEqual(len(content), info["size"])
