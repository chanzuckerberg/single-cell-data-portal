import enum

import os
import pathlib
import shutil
import tempfile

from unittest.mock import patch

import boto3
from moto import mock_s3

import backend.corpora.dataset_processing.common
import backend.corpora.dataset_processing.process_download_validate
from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
    DatasetArtifactFileType,
    UploadStatus,
    ValidationStatus,
    ProcessingStatus,
)
from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.entities.collection import Collection
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.upload import upload
from backend.corpora.common.utils.exceptions import CorporaException, MaxFileSizeExceededException
from backend.corpora.common.utils.math_utils import GB
from backend.corpora.dataset_processing.exceptions import ProcessingCancelled, ConversionFailed
from backend.corpora.dataset_processing.common import convert_file, create_artifact, get_bucket_prefix

from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestCommon(DataPortalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tmp_dir = tempfile.mkdtemp()
        cls.h5ad_filename = pathlib.Path(cls.tmp_dir, "test.h5ad")
        cls.seurat_filename = pathlib.Path(cls.tmp_dir, "test.rds")
        cls.cxg_filename = pathlib.Path(cls.tmp_dir, "test.cxg")

        cls.h5ad_filename.touch()
        cls.seurat_filename.touch()
        cls.cxg_filename.touch()

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

    def test_update_db(self):

        collection = Collection.create(self.session, visibility=CollectionVisibility.PRIVATE)
        dataset = Dataset.create(
            self.session,
            collection_id=collection.id,
        )
        dataset_id = dataset.id

        fake_env = patch.dict(os.environ, {"DATASET_ID": dataset_id, "DEPLOYMENT_STAGE": "test"})
        fake_env.start()

        backend.corpora.dataset_processing.common.update_db(
            dataset_id, metadata={"sex": ["male", "female"], "explorer_url": "https://cellxgene.com/data"}
        )
        self.session.expire(dataset)
        self.assertListEqual(Dataset.get(self.session, dataset_id).sex, ["male", "female"])

        artifact = {
            "filename": "test_filename",
            "filetype": DatasetArtifactFileType.H5AD,
            "user_submitted": True,
            "s3_uri": "s3://test_uri",
        }
        backend.corpora.dataset_processing.common.update_db(dataset_id, metadata={"artifacts": [artifact]})
        self.session.expire(dataset)
        self.assertEqual(len(Dataset.get(self.session, dataset_id).artifacts), 1)
        self.assertEqual(Dataset.get(self.session, dataset_id).artifacts[0].filename, "test_filename")
        self.assertEqual(Dataset.get(self.session, dataset_id).explorer_url, "https://cellxgene.com/data")

        backend.corpora.dataset_processing.common.update_db(
            dataset_id, processing_status={"upload_status": UploadStatus.UPLOADING, "upload_progress": 0.5}
        )
        self.session.expire(dataset)
        self.assertEqual(Dataset.get(self.session, dataset_id).processing_status.upload_status, UploadStatus.UPLOADING)
        self.assertEqual(Dataset.get(self.session, dataset_id).processing_status.upload_progress, 0.5)
        self.assertIsNone(Dataset.get(self.session, dataset_id).processing_status.validation_status)

        backend.corpora.dataset_processing.common.update_db(
            dataset_id,
            processing_status={
                "upload_status": UploadStatus.UPLOADED,
                "upload_progress": 1,
                "validation_status": ValidationStatus.VALIDATING,
            },
        )
        self.session.expire(dataset)
        self.assertEqual(Dataset.get(self.session, dataset_id).processing_status.upload_status, UploadStatus.UPLOADED)
        self.assertEqual(Dataset.get(self.session, dataset_id).processing_status.upload_progress, 1)
        self.assertEqual(
            Dataset.get(self.session, dataset_id).processing_status.validation_status, ValidationStatus.VALIDATING
        )

        fake_env.stop()

    def test_update_db__tombstoned_dataset(self):
        dataset = self.generate_dataset(self.session, tombstone=True)
        dataset_id = dataset.id

        fake_env = patch.dict(os.environ, {"DATASET_ID": dataset_id, "DEPLOYMENT_STAGE": "test"})
        fake_env.start()

        with self.assertRaises(ProcessingCancelled):
            backend.corpora.dataset_processing.common.update_db(dataset_id, metadata={"sex": ["male", "female"]})

    @patch("backend.corpora.dataset_processing.process.update_db")
    def test__create_artifact__negative(self, mock_update_db):
        artifact_bucket = "test-artifact-bucket"
        test_dataset = self.generate_dataset(
            self.session,
        )
        bucket_prefix = backend.corpora.dataset_processing.common.get_bucket_prefix(test_dataset.id)
        self.setup_s3_bucket(artifact_bucket)

        with self.subTest("file does not exist"):
            self.assertRaises(
                TypeError,
                backend.corpora.dataset_processing.common.create_artifact,
                None,
                DatasetArtifactFileType.H5AD,
                bucket_prefix,
                test_dataset.id,
                artifact_bucket,
                "h5ad_status",
            )

        with self.subTest("invalid artifact type"):

            class BadEnum(enum.Enum):
                fake = "fake"

            self.assertRaises(
                CorporaException,
                create_artifact,
                str(self.h5ad_filename),
                BadEnum.fake,
                bucket_prefix,
                test_dataset.id,
                artifact_bucket,
                "h5ad_status",
            )

        with self.subTest("dataset does not exist"):
            self.assertRaises(
                AttributeError,
                create_artifact,
                str(self.h5ad_filename),
                DatasetArtifactFileType.H5AD,
                get_bucket_prefix("1234"),
                "1234",
                artifact_bucket,
                "h5ad_status",
            )

        with self.subTest("bucket does not exist"):
            self.assertRaises(
                boto3.exceptions.S3UploadFailedError,
                create_artifact,
                str(self.h5ad_filename),
                DatasetArtifactFileType.H5AD,
                bucket_prefix,
                test_dataset.id,
                "fake-bucket",
                "h5ad_status",
            )

        # cleanup
        self.delete_s3_bucket(artifact_bucket)

    @patch("backend.corpora.dataset_processing.process.update_db")
    def test__convert_file__fail(self, mock_update_db):
        def converter(_file):
            raise RuntimeError("conversion_failed")

        with self.assertRaises(ConversionFailed):
            convert_file(converter, self.h5ad_filename, "error", "fake_id", "h5ad_status")

    def test__filesize_exceeds_limit(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection_id=collection.id)
        max_file_size_gb = CorporaConfig().upload_max_file_size_gb
        with self.assertRaises(MaxFileSizeExceededException):
            upload(
                self.session,
                collection.id,
                "download_url_placeholder",
                max_file_size_gb * GB + 1,  # Purport to exceed the max file size
                "test_user",
                dataset_id=dataset.id,
            )
        self.assertEqual(dataset.processing_status.validation_status, ValidationStatus.INVALID)
        self.assertEqual(dataset.processing_status.processing_status, ProcessingStatus.FAILURE)
        self.assertEqual(
            dataset.processing_status.validation_message,
            f"The file exceeds the maximum allowed size of {max_file_size_gb} GB",
        )
