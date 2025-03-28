import os
import unittest
from tempfile import TemporaryDirectory
from unittest import TestCase, mock

import boto3
from moto import mock_aws

from backend.common.utils.dl_sources.uri import S3URI, S3URL, CXGPublicURL, DropBoxURL, RegisteredSources, from_uri


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
        self.assertEqual("key", s3_uri.key)

    def test__validate_with_invalid_s3_uri__returns_none(self):
        s3_uri = S3URI.validate("s4://bucket/key")
        self.assertIsNone(s3_uri)

        s3_uri = S3URI.validate("s3://bucket")
        self.assertIsNone(s3_uri)

    @mock_aws
    def test__file_info__returns_file_size(self):
        s3 = boto3.client("s3")
        bucket_name = "bucket"
        key = "key"
        uri = "/".join(["s3:/", bucket_name, key])
        content = "stuff"
        s3.create_bucket(
            Bucket=bucket_name, CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]}
        )
        s3.put_object(Bucket=bucket_name, Key=key, Body=content)

        s3_uri = from_uri(uri)
        info = s3_uri.file_info()

        self.assertEqual(key, info["name"])
        self.assertEqual(len(content), info["size"])

    @mock_aws
    def test_download(self):
        s3 = boto3.client("s3")
        bucket_name = "bucket"
        key = "key/file.txt"
        uri = "/".join(["s3:/", bucket_name, key])
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
            s3_uri = from_uri(uri)
            assert isinstance(s3_uri, S3URI)
            s3_uri.download(tmpfile)
            # read the file contents
            with open(tmpfile, "r") as f:
                data = f.read()
            # assert the contents are correct
            self.assertEqual(data, content)


class TestS3URL(TestCase):

    def create_s3_url(self, bucket_name, key, content):
        s3 = boto3.client("s3")
        s3.create_bucket(
            Bucket=bucket_name, CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]}
        )
        s3.put_object(Bucket=bucket_name, Key=key, Body=content)

        return s3.generate_presigned_url("get_object", Params={"Bucket": bucket_name, "Key": key}, ExpiresIn=3600)

    @mock_aws
    def test__validate_with_valid_s3_url__ok(self):
        url = self.create_s3_url("bucket", "key/file.txt", content="stuff")
        s3_url = S3URL.validate(url)

        self.assertEqual("https", s3_url.scheme)
        self.assertEqual("bucket.s3.amazonaws.com", s3_url.netloc)
        self.assertEqual("/key/file.txt", s3_url.path)

    @mock_aws
    def test__validate_with_invalid_s3_url__returns_none(self):
        s3_url = S3URL.validate("http://somebucket.s3.amazonaws.com/key")
        self.assertIsNone(s3_url)

        s3_url = S3URL.validate("https://somebucket/key")
        self.assertIsNone(s3_url)

    @mock_aws
    def test_get_file_info(self):
        url = self.create_s3_url("bucket", "key/file.txt", content="stuff")
        s3_url = S3URL.validate(url)
        info = s3_url.file_info()
        self.assertEqual("/key/file.txt", info["name"])
        self.assertEqual(5, info["size"])


class TestCXGPubURL(TestCase):

    @mock.patch("backend.common.utils.dl_sources.uri.CorporaConfig")
    def test__validate_with_valid_url__ok(self, config_mock):
        config = mock.Mock()
        config.dataset_assets_base_url = "https://datasets.test.technology"
        config_mock.return_value = config

        url = "https://datasets.test.technology/key/file.txt"
        url = CXGPublicURL.validate(url)

        self.assertEqual("https", url.scheme)
        self.assertEqual("datasets.test.technology", url.netloc)
        self.assertEqual("/key/file.txt", url.path)

    @mock.patch("backend.common.utils.dl_sources.uri.CorporaConfig")
    def test__validate_with_invalid_url__returns_none(self, config_mock):
        config = mock.Mock()
        config.dataset_assets_base_url = "https://datasets.test.technology"
        config_mock.return_value = config

        url = CXGPublicURL.validate("http://somebucket.s3.amazonaws.com/key")
        self.assertIsNone(url)

        url = CXGPublicURL.validate("https://somebucket/key")
        self.assertIsNone(url)
