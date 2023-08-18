import os
import pathlib
import shutil
import tempfile
from unittest import TestCase

import boto3
from moto import mock_s3

from backend.common.utils.aws import delete_many_from_s3


class TestDeleteManyFromS3(TestCase):
    def setUp(self) -> None:
        super().setUp()
        self.tmp_dir = tempfile.mkdtemp()
        self.h5ad_filename = pathlib.Path(self.tmp_dir, "test.h5ad")
        self.seurat_filename = pathlib.Path(self.tmp_dir, "test.rds")
        self.cxg_filename = pathlib.Path(self.tmp_dir, "test.cxg")

        self.h5ad_filename.touch()
        self.cxg_filename.touch()
        self.seurat_filename.touch()

        # Mock S3 service if we don't have a mock api already running
        if os.getenv("BOTO_ENDPOINT_URL"):
            s3_args = {"endpoint_url": os.getenv("BOTO_ENDPOINT_URL")}
        else:
            self.mock = mock_s3()
            self.mock.start()
            s3_args = {}
        self.s3 = boto3.client("s3", config=boto3.session.Config(signature_version="s3v4"), **s3_args)

        self.uuid = "88aba04c-2d7e-4f76-8b67-071825e8bf46"
        self.bucket_name = "some-test-bucket"
        os.environ["ARTIFACT_BUCKET"] = self.bucket_name
        cxg_file = f"{self.uuid}/remixed.cxg"
        h5ad_file = f"{self.uuid}/remixed.h5ad"
        rds_file = f"{self.uuid}/remixed.rds"

        self.s3.create_bucket(
            Bucket=self.bucket_name, CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]}
        )

        self.s3.put_object(Bucket=self.bucket_name, Key=cxg_file, Body="words")
        self.s3.put_object(Bucket=self.bucket_name, Key=h5ad_file, Body="words")
        self.s3.put_object(Bucket=self.bucket_name, Key=rds_file, Body="words")

        resp = self.s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)
        assert len(resp["Contents"]) == 3
        assert resp["Contents"][0]["Key"].split("/")[0] == self.uuid

    def tearDown(self) -> None:
        self.s3.delete_bucket(Bucket=self.bucket_name)
        if not os.getenv("BOTO_ENDPOINT_URL"):
            self.mock.stop()
        shutil.rmtree(self.tmp_dir)
        super().tearDown()

    def test_delete_from_s3_deletes_all_files(self):
        resp = self.s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)
        self.assertGreater(len(resp["Contents"]), 0)
        delete_many_from_s3(self.bucket_name, self.uuid)

        resp = self.s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)
        self.assertNotIn("Contents", resp)

    def test_delete_from_s3_handles_no_files(self):
        # delete files
        delete_many_from_s3(self.bucket_name, self.uuid)
        resp = self.s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)

        self.assertNotIn("Contents", resp)
        # this should not raise any errors
        delete_many_from_s3(self.bucket_name, self.uuid)

    def test_delete_from_s3_raises_error_for_missing_dataset_id(self):
        with self.assertRaises(ValueError):
            delete_many_from_s3(self.bucket_name, "")

        response = self.s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)

        for object in response["Contents"]:
            self.s3.delete_object(Bucket=self.bucket_name, Key=object["Key"])
