import unittest

import boto3
from moto import mock_s3

from backend.corpora.common.corpora_config import CorporaConfig


class CorporaTestCaseUsingMockAWS(unittest.TestCase):
    CORPORA_TEST_CONFIG = {"bucket_name": "bogus_bucket", "region": "us-west-2"}

    def setUp(self):
        # Setup configuration
        self.corpora_config = CorporaConfig()
        self.corpora_config.set(self.__class__.CORPORA_TEST_CONFIG)

        # Setup mock AWS services
        self.s3_mock = mock_s3()
        self.s3_mock.start()

        # Corpora Bucket
        self.s3_resource = boto3.resource("s3")
        self.bucket = self.s3_resource.Bucket(self.corpora_config.bucket_name)
        self.bucket.create()

    def tearDown(self):
        super().tearDown()
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
