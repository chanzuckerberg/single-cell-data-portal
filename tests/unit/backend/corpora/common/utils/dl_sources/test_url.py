import boto3
import os
import requests
import unittest
from urllib.parse import urlparse

from backend.corpora.common.utils.dl_sources.url import from_url, RegisteredSources, S3URL, DropBoxURL
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
