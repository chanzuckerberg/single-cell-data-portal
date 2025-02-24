from unittest.mock import Mock, patch

import pytest

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.common.entities import (
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetUploadStatus,
    DatasetValidationStatus,
)
from backend.layers.common.ingestion_manifest import IngestionManifest
from backend.layers.processing.process import ProcessMain
from backend.layers.processing.process_validate_h5ad import ProcessValidateH5AD
from tests.unit.processing.base_processing_test import BaseProcessingTest


class TestProcessDownload(BaseProcessingTest):
    @patch("backend.common.utils.dl_sources.uri.S3Provider")
    @patch("backend.common.utils.dl_sources.uri.S3URI.file_info", return_value={"size": 100, "name": "fake_name"})
    def test_download_from_s3_uri(self, *arg):
        """
        Call process download using an s3 uri
        """

        s3_uri = "s3://fake_bucket_name/fake_key/fake_file.h5ad"
        pdv = ProcessMain(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pdv.download_from_s3 = Mock()

        assert pdv.download_from_source_uri(s3_uri, "fake_local_path") == "fake_local_path"

    @patch("backend.common.utils.dl_sources.uri.downloader")
    @patch("backend.common.utils.dl_sources.uri.DropBoxURL.file_info", return_value={"size": 100, "name": "fake_name"})
    def test_download_from_dropbox_uri(self, *arg):
        """
        Call process download using a dropbox uri
        """

        dropbox_uri = "https://www.dropbox.com/s/fake_location/test.h5ad?dl=1"
        pdv = ProcessMain(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pdv.download = Mock()

        assert pdv.download_from_source_uri(dropbox_uri, "fake_local_path") == "fake_local_path"

    def test_download_unknown_uri(self):
        """
        Call process download using unknown
        """

        uri = "fake://fake_bucket_name/fake_key/fake_file.h5ad"
        pdv = ProcessMain(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pdv.download_from_s3 = Mock()
        with pytest.raises(ValueError, match=f"Malformed source URI: {uri}"):
            pdv.download_from_source_uri(uri, "fake_local_path")


class TestProcessValidateH5AD(BaseProcessingTest):
    def test_process_download_validate_success(self):
        """
        ProcessDownloadValidate should:
        1. Download the h5ad artifact
        2. Set DatasetStatusKey.H5AD DatasetValidationStatus.VALIDATING
        3. Validate the h5ad
        4. Set DatasetStatusKey.H5AD DatasetValidationStatus.VALID
        5. Set the DatasetStatusKey.RDS DatasetConversionStatus.SKIPPED accordingly
        6. upload the original file to S3

        """
        dropbox_uri = "https://www.dropbox.com/s/fake_location/test.h5ad?dl=0"

        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )
        # This is where we're at when we start the SFN

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertIsNone(status.validation_status)
        self.assertEqual(status.processing_status, DatasetProcessingStatus.INITIALIZED)
        self.assertEqual(status.upload_status, DatasetUploadStatus.WAITING)

        pdv = ProcessValidateH5AD(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pdv.download_from_source_uri = Mock(return_value=CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME)
        pdv.process(dataset_version_id, IngestionManifest(anndata=dropbox_uri), "fake_bucket_name")
        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.rds_status, DatasetConversionStatus.SKIPPED)
        self.assertEqual(status.h5ad_status, DatasetConversionStatus.CONVERTING)

        # Verify that the original (raw.h5ad) file is there
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/raw.h5ad"))
        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
        self.assertEqual(1, len(artifacts))

    def test_process_validate_fail(self):
        """
        If the validation is not successful, the processing pipeline should:
        1. Set the processing status to INVALID
        2. Set a validation message accordingly
        """
        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        manifest = IngestionManifest(anndata=dropbox_uri)
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )

        # Set a mock failure for the schema validator
        self.schema_validator.validate_anndata = Mock(
            return_value=(False, ["Validation error 1", "Validation error 2"], True)
        )
        pm = ProcessMain(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pm.download_from_source_uri = Mock(return_value=CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME)
        for step_name in ["validate_anndata"]:
            pm.process(
                collection.version_id,
                dataset_version_id,
                step_name,
                manifest,
                "fake_bucket_name",
                "fake_datasets_bucket",
                "fake_cxg_bucket",
            )

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.validation_status, DatasetValidationStatus.INVALID)
        self.assertEqual(status.validation_anndata_message, "Validation error 1\nValidation error 2")
