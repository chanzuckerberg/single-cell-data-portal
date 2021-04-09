import boto3
import os
import pathlib
import requests
import shutil
import tempfile
from moto import mock_s3
from unittest.mock import patch

from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
)
from backend.corpora.dataset_processing import process
from backend.corpora.dataset_processing.process import make_cxg, make_seurat, make_loom
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestDatasetProcessing(DataPortalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tmp_dir = tempfile.mkdtemp()
        cls.h5ad_filename = pathlib.Path(cls.tmp_dir, "test.h5ad")
        cls.seurat_filename = pathlib.Path(cls.tmp_dir, "test.rds")
        cls.loom_filename = pathlib.Path(cls.tmp_dir, "test.loom")
        cls.cxg_filename = pathlib.Path(cls.tmp_dir, "test.cxg")
        cls.real_h5ad_filename = pathlib.Path(cls.tmp_dir, "real.h5ad")
        cls.h5ad_filename.touch()
        cls.seurat_filename.touch()
        cls.loom_filename.touch()
        cls.cxg_filename.touch()

        presigned_url = cls.get_presigned_url(
            "a5da03d7-f44b-4bb8-81a4-bc80960d6223", "04be802c-ffc6-4a7f-a4c5-f2812f1a707f"
        )
        cls.download(presigned_url, cls.real_h5ad_filename)

    @staticmethod
    def get_presigned_url(dataset_id, asset_id):
        response = requests.post(
            f"https://api.cellxgene.staging.single-cell.czi.technology/dp/v1/datasets/" f"{dataset_id}/asset/{asset_id}"
        )
        response.raise_for_status()
        return response.json()["presigned_url"]

    @staticmethod
    def download(url, local_filename):
        with requests.get(url, stream=True) as resp:
            resp.raise_for_status()
            with open(local_filename, "wb") as fp:
                for chunk in resp.iter_content(chunk_size=None):
                    fp.write(chunk)

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        shutil.rmtree(cls.tmp_dir)

    def setup_s3_bucket(self, bucket_name):
        # Mock S3 service if we don't have a mock api already running
        if os.getenv("BOTO_ENDPOINT_URL"):
            s3_args = {"endpoint_url": os.getenv("BOTO_ENDPOINT_URL")}
        else:
            s3_mock = mock_s3()
            s3_mock.start()
            s3_args = {}
            self.addCleanup(s3_mock.stop)
        s3 = boto3.client("s3", config=boto3.session.Config(signature_version="s3v4"), **s3_args)
        self.s3_resource = boto3.resource("s3", config=boto3.session.Config(signature_version="s3v4"), **s3_args)
        try:
            s3.create_bucket(
                Bucket=bucket_name, CreateBucketConfiguration={"LocationConstraint": os.environ["AWS_DEFAULT_REGION"]}
            )
        except self.s3_resource.meta.client.exceptions.BucketAlreadyExists:
            pass
        return s3

    def delete_s3_bucket(self, bucket_name):
        bucket = self.s3_resource.Bucket(bucket_name)
        if bucket.creation_date is not None:
            bucket.objects.all().delete()
            bucket.delete()

    def test_make_cxg(self):
        make_cxg(str(self.real_h5ad_filename))

    def test_make_seurat(self):
        make_seurat(str(self.real_h5ad_filename))

    def test_make_loom(self):
        make_loom(str(self.real_h5ad_filename))

    @patch("backend.corpora.dataset_processing.process.get_download_url_from_shared_link")
    def test_main(self, mock_get_download_url_from_shared_link):
        url = self.get_presigned_url("a5da03d7-f44b-4bb8-81a4-bc80960d6223", "04be802c-ffc6-4a7f-a4c5-f2812f1a707f")
        mock_get_download_url_from_shared_link.return_value = url
        dataset = self.generate_dataset(
            self.session, collection_id="test_collection_id", collection_visibility=CollectionVisibility.PUBLIC.name
        )
        process.process(dataset.id, url, "s3://cellxgene", "s3://artifacts")
