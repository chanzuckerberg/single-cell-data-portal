import os
import unittest

import boto3
from moto import mock_s3

from backend.corpora.common.corpora_config import CorporaConfig


class CorporaTestCaseUsingMockAWS(unittest.TestCase):
    CORPORA_TEST_CONFIG = {"bucket_name": "bogus-bucket"}

    def setUp(self):
        # Setup configuration
        self.corpora_config = CorporaConfig()
        self.corpora_config.set(self.__class__.CORPORA_TEST_CONFIG)

        # Mock S3 service if we don't have a mock api already running
        if not os.getenv("BOTO_ENDPOINT_URL"):
            self.s3_mock = mock_s3()
            self.s3_mock.start()

        # Corpora Bucket
        self.s3_resource = boto3.resource(
            "s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"), config=boto3.session.Config(signature_version="s3v4")
        )
        self.bucket = self.s3_resource.Bucket(self.corpora_config.bucket_name)
        self.bucket.create()

    def tearDown(self):
        super().tearDown()
        if not os.getenv("BOTO_ENDPOINT_URL"):
            self.s3_mock.stop()

    def create_s3_object(
        self, object_key, bucket_name=None, content_type="application/octet-stream", content="file_content"
    ):
        bucket_name = bucket_name or self.bucket.name

        if not self.s3_resource.Bucket(bucket_name) in self.s3_resource.buckets.all():
            self.s3_resource.Bucket(bucket_name).create()

        s3object = self.s3_resource.Bucket(bucket_name).Object(object_key)
        s3object.put(Body=content, ContentType=content_type)
        return s3object
