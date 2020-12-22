import os
import pathlib
import shutil
from unittest import TestCase
import tempfile
import boto3
from moto import mock_s3

from backend.corpora.lambdas.upload_failures.upload import delete_many_from_s3

s3 = boto3.client("s3")


class TestUploadFailureHandling(TestCase):
    def setUp(self) -> None:
        self.tmp_dir = tempfile.mkdtemp()
        self.h5ad_filename = pathlib.Path(self.tmp_dir, "test.h5ad")
        self.seurat_filename = pathlib.Path(self.tmp_dir, "test.rds")
        self.loom_filename = pathlib.Path(self.tmp_dir, "test.loom")
        self.cxg_filename = pathlib.Path(self.tmp_dir, "test.cxg")

        self.h5ad_filename.touch()
        self.cxg_filename.touch()
        self.seurat_filename.touch()
        self.loom_filename.touch()

        self.mock = mock_s3()
        self.mock.start()
        self.uuid = "88aba04c-2d7e-4f76-8b67-071825e8bf46"
        self.bucket_name = "some_Bucket"
        os.environ["ARTIFACT_BUCKET"] = self.bucket_name
        cxg_file = f"{self.uuid}/remixed.cxg"
        h5ad_file = f"{self.uuid}/remixed.h5ad"
        rds_file = f"{self.uuid}/remixed.rds"
        loom_file = f"{self.uuid}/remixed.loom"
        s3.create_bucket(
            Bucket=self.bucket_name, CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]}
        )

        s3.put_object(Bucket=self.bucket_name, Key=cxg_file, Body="words")
        s3.put_object(Bucket=self.bucket_name, Key=h5ad_file, Body="words")
        s3.put_object(Bucket=self.bucket_name, Key=rds_file, Body="words")
        s3.put_object(Bucket=self.bucket_name, Key=loom_file, Body="words")

        resp = s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)
        assert len(resp["Contents"]) == 4
        assert resp["Contents"][0]["Key"].split("/")[0] == self.uuid

    def tearDown(self) -> None:
        self.mock.stop()
        shutil.rmtree(self.tmp_dir)

    def test_delete_from_s3_deletes_all_files(self):
        resp = s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)
        self.assertGreater(len(resp["Contents"]), 0)
        delete_many_from_s3(self.bucket_name, self.uuid)

        resp = s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)
        self.assertNotIn("Contents", resp)

    def test_delete_from_s3_handles_no_files(self):
        # delete files
        delete_many_from_s3(self.bucket_name, self.uuid)
        resp = s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)

        self.assertNotIn("Contents", resp)
        # this should not raise any errors
        delete_many_from_s3(self.bucket_name, self.uuid)
