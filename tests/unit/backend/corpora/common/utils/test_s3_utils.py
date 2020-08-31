import unittest

import boto3
from moto import mock_s3

from backend.corpora.common.utils.s3_utils import generate_file_url, head_file


class TestS3Utils(unittest.TestCase):
    def setUp(self):
        self.s3_mock = mock_s3()
        self.s3_mock.start()

    def tearDown(self):
        self.s3_mock.stop()

    def test_generate_file_url(self):
        file_prefix = "test_prefix"
        bucket_name = "mock-bucket"

        with self.subTest("Generate URL OK"):
            url = generate_file_url(bucket_name, file_prefix)
            self.assertIn(file_prefix, url)
            self.assertIn("Expires=", url)

    def test_head_file(self):
        file_name = "fake_file.txt"
        body = "hello world!"
        bucket = "corpora-data-test"
        region = "us-west-2"
        conn = boto3.resource("s3", region_name=region)
        conn.create_bucket(Bucket=bucket)
        s3 = boto3.client("s3", region_name=region)
        s3.put_object(Bucket=bucket, Key=file_name, Body=body)

        response = head_file(bucket, file_name)
        self.assertEqual(len(body), response["ContentLength"])

    def test_head_file__not_found(self):
        file_name = "fake_file.txt"
        bucket = "corpora-data-test"
        region = "us-west-2"
        conn = boto3.resource("s3", region_name=region)
        conn.create_bucket(Bucket=bucket)
        s3 = boto3.client("s3", region_name=region)
        self.assertIs(None, head_file(bucket, file_name, s3=s3))
