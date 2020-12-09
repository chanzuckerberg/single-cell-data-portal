import os
import pathlib
import shutil
from unittest import TestCase
import tempfile
import boto3
from moto import mock_s3

from backend.corpora.lambdas.upload_failures.upload import delete_many_from_s3

s3 = boto3.client("s3")


class TestUploadFailure(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp_dir = tempfile.mkdtemp()
        cls.h5ad_filename = pathlib.Path(cls.tmp_dir, "test.h5ad")
        cls.seurat_filename = pathlib.Path(cls.tmp_dir, "test.rds")
        cls.loom_filename = pathlib.Path(cls.tmp_dir, "test.loom")
        cls.cxg_filename = pathlib.Path(cls.tmp_dir, "test.cxg")

        cls.h5ad_filename.touch()
        cls.cxg_filename.touch()
        cls.seurat_filename.touch()
        cls.loom_filename.touch()

        cls.mock = mock_s3()
        cls.mock.start()
        cls.uuid = "88aba04c-2d7e-4f76-8b67-071825e8bf46"
        cls.bucket_name = "some_Bucket"
        os.environ["ARTIFACT_BUCKET"] = cls.bucket_name
        cxg_file = f'{cls.uuid}/remixed.cxg'
        h5ad_file = f'{cls.uuid}/remixed.h5ad'
        rds_file = f'{cls.uuid}/remixed.rds'
        loom_file = f'{cls.uuid}/remixed.loom'
        s3.create_bucket(Bucket=cls.bucket_name)

        s3.put_object(Bucket=cls.bucket_name, Key=cxg_file, Body="words")
        s3.put_object(Bucket=cls.bucket_name, Key=h5ad_file, Body="words")
        s3.put_object(Bucket=cls.bucket_name, Key=rds_file, Body="words")
        s3.put_object(Bucket=cls.bucket_name, Key=loom_file, Body="words")

        resp = s3.list_objects_v2(Bucket=cls.bucket_name, Prefix=cls.uuid)
        assert len(resp['Contents']) == 4
        assert resp['Contents'][0]['Key'].split('/')[0] == cls.uuid

    @classmethod
    def tearDownClass(cls) -> None:
        cls.mock.stop()
        shutil.rmtree(cls.tmp_dir)

    def test_delete_from_s3_deletes_all_files(self):
        resp = s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)
        self.assertGreater(len(resp['Contents']), 0)
        delete_many_from_s3(self.bucket_name, self.uuid)

        resp = s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)
        self.assertNotIn('Contents', resp)

    def test_delete_from_s3_handles_no_files(self):
        # delete files
        # delete_many_from_s3(self.bucket_name, self.uuid)
        resp = s3.list_objects_v2(Bucket=self.bucket_name, Prefix=self.uuid)

        self.assertNotIn('Contents', resp)
        # this should not raise any errors
        delete_many_from_s3(self.bucket_name, self.uuid)


