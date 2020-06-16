import unittest

import boto3
from moto import mock_iam, mock_s3

from backend.corpora.common.corpora_config import CorporaConfig


class CorporaTestCaseUsingMockAWS(unittest.TestCase):
    CORPORA_TEST_CONFIG = {
        "bucket_name": "bogus_bucket"
    }

    def setUp(self):
        # Setup configuration
        self.corpora_config = CorporaConfig()
        self.corpora_config.set(self.__class__.CORPORA_TEST_CONFIG)

        # Setup mock AWS services
        self.s3_mock = mock_s3()
        self.s3_mock.start()
        self.iam_mock = mock_iam()
        self.iam_mock.start()

        # Corpora Bucket
        self.s3_resource = boto3.resource("s3")
        self.bucket = self.s3_resource.Bucket(self.corpora_config.bucket_name)
        self.bucket.create()

    def tearDown(self):
        super().tearDown()
        self.s3_mock.stop()
        self.iam_mock.stop()

    def create_s3_object(self, object_key, bucket_name=None, content_type="application/octet-stream",
                         content="file_content"):
        bucket_name = bucket_name or self.bucket.name
        s3object = self.s3_resource.Bucket(bucket_name).Object(object_key)
        s3object.put(Body=content, ContentType=content_type)
        return s3object
