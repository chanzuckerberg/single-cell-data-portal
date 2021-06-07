import random
import tempfile
import os

import botocore
import boto3
import typing
from boto.s3.bucket import Bucket
from moto import mock_s3

from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.corpora_orm import DatasetArtifactType, DatasetArtifactFileType
from backend.corpora.common.entities import DatasetAsset, Dataset
from backend.corpora.common.entities.dataset import get_cxg_bucket_path
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

    def generate_artifact(
        self, session, dataset_id, artifact_type=DatasetArtifactFileType.H5AD, file_name="data", upload=False
    ) -> DatasetAsset:
        file_name = f"{file_name}.{artifact_type.name}"
        if upload:
            with tempfile.TemporaryDirectory() as temp_path:
                temp_file = f"{temp_path}/{file_name}"
                content = "".join(random.choices("abcdef", k=16))
                with open(temp_file, "w") as fp:
                    fp.write(content)
                s3_uri = DatasetAsset.upload(temp_file, dataset_id, self.bucket.name)
        else:
            s3_uri = DatasetAsset.make_s3_uri(self.bucket.name, dataset_id, file_name)
        return DatasetAsset.create(
            session, dataset_id, file_name, artifact_type, DatasetArtifactType.REMIX, False, s3_uri
        )

    def create_explorer_s3_object(self, session, dataset_id, upload=False):
        if upload:
            file_name = f"{dataset_id}.cxg/"
            content = "".join(random.choices("abcdef", k=16))
            self.cellxgene_bucket.Object(file_name).put(Body=content, ContentType="application/octet-stream")
        dataset = Dataset.get(session, dataset_id)
        explorer_url = f"http://bogus.url/d/{dataset_id}.cxg/"
        # explorer_s3_uri = f"s3://{self.cellxgene_bucket}/{dataset_id}.cxg/" TODO use line once cxg in DatasetArtifact
        dataset.update(explorer_url=explorer_url)

    def generate_dataset_with_s3_resources(self, session, artifacts=True, explorer_s3_object=True, **params) -> Dataset:
        dataset = self.generate_dataset(session, **params)
        if artifacts:
            for ext in DatasetArtifactFileType:
                self.generate_artifact(session, dataset.id, ext, upload=True)
        if explorer_s3_object:
            self.create_explorer_s3_object(session, dataset.id, upload=True)
        return dataset

    def get_s3_object_paths_from_dataset(self, dataset: Dataset) -> typing.List[typing.Tuple[Bucket, str]]:
        s3_objects = [(self.bucket, DatasetAsset(art).get_bucket_path()) for art in dataset.artifacts]
        if dataset.explorer_url:
            s3_objects + [(self.cellxgene_bucket, f"{get_cxg_bucket_path(dataset.explorer_url)}.cxg/")]
        return s3_objects

    def assertS3FileExists(self, bucket: Bucket, file_name: str):
        self.assertGreater(bucket.Object(file_name).content_length, 1)

    def assertS3FileDoesNotExist(self, bucket: Bucket, file_name: str, msg: str = None):
        msg = msg if msg else f"s3://{bucket.name}/{file_name} found."
        with self.assertRaises(botocore.exceptions.ClientError, msg=msg):
            bucket.Object(file_name).content_length
