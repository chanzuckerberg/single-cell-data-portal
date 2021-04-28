import os

import boto3
from moto import mock_s3

from backend.corpora.common.corpora_config import CorporaConfig
from tests.unit.backend.fixtures import config
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class CorporaTestCaseUsingMockAWS(DataPortalTestCase):
    CORPORA_TEST_CONFIG = config.CORPORA_TEST_CONFIG

    def setUp(self):
        super().setUp()
        # Setup configuration
        self.corpora_config = CorporaConfig()
        self.corpora_config.set(config.CORPORA_TEST_CONFIG)

        # Mock S3 service if we don't have a mock api already running
        if os.getenv("BOTO_ENDPOINT_URL"):
            s3_args = {"endpoint_url": os.getenv("BOTO_ENDPOINT_URL")}
        else:
            self.s3_mock = mock_s3()
            self.s3_mock.start()
            s3_args = {}

        self.s3_resource = boto3.resource("s3", config=boto3.session.Config(signature_version="s3v4"), **s3_args)
        # Corpora Bucket
        self.bucket = self.s3_resource.Bucket(self.corpora_config.bucket_name)
        try:
            self.bucket.create(CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]})
        except self.s3_resource.meta.client.exceptions.BucketAlreadyExists:
            pass

        # Cellxgene Bucket
        self.cellxgene_bucket = self.s3_resource.Bucket(
            os.getenv("CELLXGENE_BUCKET", f"hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}")
        )
        try:
            self.cellxgene_bucket.create(
                CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]}
            )
        except self.s3_resource.meta.client.exceptions.BucketAlreadyExists:
            pass

    def tearDown(self):
        super().tearDown()
        self.bucket.objects.all().delete()
        self.cellxgene_bucket.objects.all().delete()
        self.cellxgene_bucket.delete()
        self.bucket.delete()
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
