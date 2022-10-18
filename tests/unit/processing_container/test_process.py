import os
import shutil
from unittest.mock import patch

import requests

from backend.common.corpora_orm import (
    ConversionStatus,
    UploadStatus,
    ValidationStatus,
    DatasetArtifactFileType,
)
from backend.common.entities import Dataset
from backend.dataset_processing.exceptions import (
    ProcessingFailed,
    ValidationFailed,
    ConversionFailed,
)
from backend.dataset_processing.process import main
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from tests.unit.backend.fixtures.environment_setup import EnvironmentSetup


class TestDatasetProcessing(CorporaTestCaseUsingMockAWS):
    """
    This test class is intended to exercise the upload pipeline processing, running within
    a "corpora-upload" (upload processing) Docker container. As the Docker container is configured
    with all necessary software dependencies, manually running this test outside of Docker
    (e.g. within an IDE when making changes to the test itself) will require you to:
    * set env vars: CORPORA_LOCAL_DEV=1;BOTO_ENDPOINT_URL=http://localhost:4566
    * Locally install latest cellxgene-schema (via pip)
    * Install all the same Python dependencies at Dockerfile.processing_image:23 and Dockerfile.processing_base:19
      (R software) to ensure Seurat artifact can be generated. (These will only matter if we add test
      assertions to check these artifact's creation.)
    """

    @staticmethod
    def fixture_file_path(relative_filename):
        return os.path.abspath(os.path.join(os.path.dirname(__file__), relative_filename))

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.h5ad_raw = cls.fixture_file_path("fixtures/3_0_0_valid.h5ad")

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
        cls.clean_generated_files()

    @classmethod
    def clean_generated_files(cls):
        if os.path.exists("local.cxg"):
            shutil.rmtree("local.cxg")
        for f in ["local.h5ad", "local.rds"]:
            if os.path.exists(f):
                os.remove(f)

    @patch("backend.dataset_processing.process_download_validate.download_from_source_uri")
    @patch("backend.dataset_processing.remaster_cxg.process")  # TODO: provide test data to properly test this.
    def test_main(self, mock_remaster_cxg_process, mock_download_from_source_uri):
        """
        Tests full pipeline for processing an uploaded H5AD file, including database updates
        generation and upload of all artifacts to S3 (localstack), but excluding the Dropbox download
        functionality.  Dropbox I/O is mocked to prevent dependency on remote services (non-Dockerized).
        """
        mock_download_from_source_uri.return_value = self.h5ad_raw
        dataset = self.generate_dataset(self.session, collection_id="test_collection_id")
        test_environment = {
            "STEP_NAME": "cxg",
            "DROPBOX_URL": "https://www.dropbox.com/IGNORED",
            "ARTIFACT_BUCKET": self.corpora_config.bucket_name,
            "CELLXGENE_BUCKET": self.corpora_config.bucket_name,
            "DATASET_ID": dataset.id,
            "DEPLOYMENT_STAGE": "test",
        }

        test_environment["STEP_NAME"] = "download-validate"
        with EnvironmentSetup(test_environment):
            main()
            artifacts = [artifact.filetype == DatasetArtifactFileType.H5AD for artifact in dataset.artifacts]
            self.assertTrue(all(artifacts), "Assert all H5AD")
            self.assertEqual(2, len(artifacts), "Assert 2 H5ADs")

        test_environment["STEP_NAME"] = "cxg"
        with EnvironmentSetup(test_environment):
            main()

        test_environment["STEP_NAME"] = "seurat"
        with EnvironmentSetup(test_environment):
            main()

        test_environment["STEP_NAME"] = "cxg_remaster"
        with EnvironmentSetup(test_environment):
            main()
            mock_remaster_cxg_process.assert_called()

        # TODO: add assertions. See https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data-portal/1449 # noqa: E501
        # 1. H5AD has annotation labels added and uploaded to S3
        # 2. cxg, rds uploaded to s3
        # 3. databases metadata updated and showing successful status

    @patch("backend.dataset_processing.process_cxg.process")
    @patch("backend.dataset_processing.process_seurat.process")
    @patch("backend.dataset_processing.process_download_validate.process")
    def test_main__Exceptions(self, mock_process_download_validate, mock_process_seurat, mock_process_cxg):
        test_environment = {
            "DROPBOX_URL": "https://www.dropbox.com/IGNORED",
            "ARTIFACT_BUCKET": self.corpora_config.bucket_name,
            "CELLXGENE_BUCKET": self.corpora_config.bucket_name,
            "DEPLOYMENT_STAGE": "test",
        }

        test_environment["STEP_NAME"] = "download-validate"
        with self.subTest("download-validate step raises unknown exception"):
            mock_process_download_validate.side_effect = Exception("unknown failure")
            dataset = self.generate_dataset(self.session, collection_id="test_collection_id")
            test_environment["DATASET_ID"] = dataset.id
            with EnvironmentSetup(test_environment):
                main()
                expected_dataset = Dataset.get(self.session, test_environment["DATASET_ID"]).processing_status
                self.assertEqual(expected_dataset.upload_status, UploadStatus.FAILED)
                self.assertEqual(expected_dataset.upload_message, "unknown failure")

        with self.subTest("download-validate step raises ValidationFailed"):
            mock_process_download_validate.side_effect = ValidationFailed(
                {"validation_status": ValidationStatus.INVALID, "validation_message": "Schema version is not supported"}
            )
            dataset = self.generate_dataset(self.session, collection_id="test_collection_id")
            test_environment["DATASET_ID"] = dataset.id
            with EnvironmentSetup(test_environment):
                main()
                expected_dataset = Dataset.get(self.session, test_environment["DATASET_ID"]).processing_status
                self.assertEqual(expected_dataset.validation_status, ValidationStatus.INVALID)
                self.assertEqual(expected_dataset.validation_message, "Schema version is not supported")

        with self.subTest("download-validate step raises ProcessingFailed"):
            mock_process_download_validate.side_effect = ProcessingFailed(
                {"upload_status": UploadStatus.FAILED, "upload_message": "Insufficient disk space"}
            )
            dataset = self.generate_dataset(self.session, collection_id="test_collection_id")
            test_environment["DATASET_ID"] = dataset.id
            with EnvironmentSetup(test_environment):
                main()
                expected_dataset = Dataset.get(self.session, test_environment["DATASET_ID"]).processing_status
                self.assertEqual(expected_dataset.upload_status, UploadStatus.FAILED)
                self.assertEqual(expected_dataset.upload_message, "Insufficient disk space")

        with self.subTest("seurat step raises unknown exception"):
            mock_process_seurat.side_effect = Exception("unknown failure")
            dataset = self.generate_dataset(self.session, collection_id="test_collection_id")
            test_environment["DATASET_ID"] = dataset.id
            test_environment["STEP_NAME"] = "seurat"
            with EnvironmentSetup(test_environment):
                main()
                expected_dataset = Dataset.get(self.session, test_environment["DATASET_ID"]).processing_status
                self.assertEqual(expected_dataset.rds_status, ConversionStatus.FAILED)

        with self.subTest("cxg step raises unknown exception"):
            mock_process_cxg.side_effect = Exception("unknown failure")
            dataset = self.generate_dataset(self.session, collection_id="test_collection_id")
            test_environment["DATASET_ID"] = dataset.id
            test_environment["STEP_NAME"] = "cxg"
            with EnvironmentSetup(test_environment):
                main()
                expected_dataset = Dataset.get(self.session, test_environment["DATASET_ID"]).processing_status
                self.assertEqual(expected_dataset.cxg_status, ConversionStatus.FAILED)

        with self.subTest("conversion step raises ConversionFailed"):
            mock_process_cxg.side_effect = ConversionFailed({"cxg_status": ConversionStatus.FAILED})
            dataset = self.generate_dataset(self.session, collection_id="test_collection_id")
            test_environment["DATASET_ID"] = dataset.id
            test_environment["STEP_NAME"] = "cxg"
            with EnvironmentSetup(test_environment):
                main()
                expected_dataset = Dataset.get(self.session, test_environment["DATASET_ID"]).processing_status
                self.assertEqual(expected_dataset.cxg_status, ConversionStatus.FAILED)
