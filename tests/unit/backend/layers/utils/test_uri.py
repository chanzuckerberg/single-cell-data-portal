import os
import unittest
from tempfile import TemporaryDirectory
from unittest import TestCase

import boto3
from moto import mock_s3

from backend.common.utils.dl_sources.uri import S3URI, S3URL, DropBoxURL, RegisteredSources, from_uri


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

    def test_from_uri(self):
        RegisteredSources.add(S3URL)
        self.assertIsInstance(from_uri("https://somebucket.s3.amazonaws.com"), S3URL)

        RegisteredSources.add(DropBoxURL)
        self.assertIsInstance(from_uri("https://www.dropbox.com/s/abcd1234/abcd1234?dl=0"), DropBoxURL)
        self.assertIsNone(from_uri("http://example.com"))

        RegisteredSources.add(S3URI)
        self.assertIsInstance(from_uri("s3://bucket/key"), S3URI)


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

        s3_uri = from_uri(f"s3://{bucket_name}{key}")
        info = s3_uri.file_info()

        self.assertEqual(key, info["name"])
        self.assertEqual(len(content), info["size"])

    @mock_s3
    def test_download(self):
        s3 = boto3.client("s3")
        bucket_name = "bucket"
        key = "/key"
        content = "stuff"
        s3.create_bucket(
            Bucket=bucket_name, CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]}
        )
        s3.put_object(Bucket=bucket_name, Key=key, Body=content)

        # add a context manager for a temporary directory
        with TemporaryDirectory() as tmpdir:
            # create a temporary file
            tmpfile = os.path.join(tmpdir, "test.txt")
            # download the file to the temporary directory
            s3_uri = from_uri(f"s3://{bucket_name}{key}")
            s3_uri.download(tmpfile)
            # read the file contents
            with open(tmpfile, "r") as f:
                data = f.read()
            # assert the contents are correct
            self.assertEqual(data, content)
