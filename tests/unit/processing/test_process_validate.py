import tempfile
from unittest.mock import Mock, patch

import anndata
import pytest

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.common.entities import (
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersionId,
    Link,
)
from backend.layers.common.ingestion_manifest import to_manifest
from backend.layers.processing.process import ProcessMain
from backend.layers.processing.process_validate import ProcessValidate
from tests.unit.processing.base_processing_test import BaseProcessingTest


class TestProcessDownload(BaseProcessingTest):
    @patch("backend.common.utils.dl_sources.uri.S3Provider")
    @patch("backend.common.utils.dl_sources.uri.S3URI.file_info", return_value={"size": 100, "name": "fake_name"})
    def test_download_from_s3_uri(self, *arg):
        """
        Call process download using an s3 uri
        """

        s3_uri = "s3://fake_bucket_name/fake_key/fake_file.h5ad"
        pdv = ProcessValidate(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pdv.download_from_s3 = Mock()

        assert pdv.download_from_source_uri(s3_uri, "fake_local_path") == "fake_local_path"

    @patch("backend.common.utils.dl_sources.uri.downloader")
    @patch("backend.common.utils.dl_sources.uri.DropBoxURL.file_info", return_value={"size": 100, "name": "fake_name"})
    def test_download_from_dropbox_uri(self, *arg):
        """
        Call process download using a dropbox uri
        """

        dropbox_uri = "https://www.dropbox.com/s/fake_location/test.h5ad?dl=1"
        pdv = ProcessValidate(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pdv.download = Mock()

        assert pdv.download_from_source_uri(dropbox_uri, "fake_local_path") == "fake_local_path"

    def test_download_unknown_uri(self):
        """
        Call process download using unknown
        """

        uri = "fake://fake_bucket_name/fake_key/fake_file.h5ad"
        pdv = ProcessValidate(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        pdv.download_from_s3 = Mock()
        with pytest.raises(ValueError, match=f"Malformed source URI: {uri}"):
            pdv.download_from_source_uri(uri, "fake_local_path")


class ProcessingTest(BaseProcessingTest):
    def test_process_download_validate_success(self):
        """
        ProcessDownloadValidate should:
        1. Download the h5ad artifact
        2. set validation status to VALID
        3. Set upload status to UPLOADED
        4. set h5ad status to UPLOADED
        5. upload the original file to S3
        6. upload the labeled file to S3
        """
        dropbox_uri = "https://www.dropbox.com/s/fake_location/test.h5ad?dl=0"

        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )
        # This is where we're at when we start the SFN

        status = self.business_logic.get_dataset_status(dataset_version_id)
        # self.assertEqual(status.validation_status, DatasetValidationStatus.NA)
        self.assertIsNone(status.validation_status)
        self.assertEqual(status.processing_status, DatasetProcessingStatus.INITIALIZED)
        self.assertEqual(status.upload_status, DatasetUploadStatus.WAITING)

        with patch("backend.layers.processing.process_validate.ProcessValidate.extract_metadata"):
            pdv = ProcessValidate(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
            pdv.download_from_source_uri = Mock(return_value=CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME)
            pdv.populate_dataset_citation = Mock()
            pdv.process(
                collection.version_id, dataset_version_id, "fake_uri", "fake_bucket_name", "fake_datasets_bucket"
            )
            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.validation_status, DatasetValidationStatus.VALID)
            self.assertEqual(status.h5ad_status, DatasetConversionStatus.UPLOADED)
            pdv.populate_dataset_citation.assert_called_once_with(
                collection.version_id, dataset_version_id, CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME
            )

            # Verify that both the original (raw.h5ad) and the labeled (local.h5ad) files are there
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/raw.h5ad"))
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"))
            # Verify that the labeled file is uploaded to the datasets bucket
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}.h5ad"))

            artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
            self.assertEqual(2, len(artifacts))

    def test_populate_dataset_citation__with_publication_doi(self):
        mock_adata = anndata.AnnData(X=None, obs=None, obsm=None, uns={}, var=None)
        self.crossref_provider.fetch_metadata = Mock(return_value=({}, "12.2345", 17169328.664))
        collection = self.generate_unpublished_collection(
            links=[Link(name=None, type="DOI", uri="https://doi.org/12.2345")]
        )
        with tempfile.NamedTemporaryFile(suffix=".h5ad") as f:
            mock_adata.write_h5ad(f.name)
            pdv = ProcessValidate(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
            dataset_version_id = DatasetVersionId()
            pdv.populate_dataset_citation(collection.version_id, dataset_version_id, f.name)
            citation_str = (
                f"Publication: https://doi.org/12.2345 "
                f"Dataset Version: http://domain/{dataset_version_id}.h5ad curated and distributed by "
                f"CZ CELLxGENE Discover in Collection: https://domain/collections/{collection.collection_id}"
            )
            adata = anndata.read_h5ad(f.name)
            self.assertEqual(adata.uns["citation"], citation_str)

    def test_populate_dataset_citation__no_publication_doi(self):
        mock_adata = anndata.AnnData(X=None, obs=None, obsm=None, uns={}, var=None)
        collection = self.generate_unpublished_collection()
        with tempfile.NamedTemporaryFile(suffix=".h5ad") as f:
            mock_adata.write_h5ad(f.name)
            pdv = ProcessValidate(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
            dataset_version_id = DatasetVersionId()
            pdv.populate_dataset_citation(collection.version_id, dataset_version_id, f.name)
            citation_str = (
                f"Dataset Version: http://domain/{dataset_version_id}.h5ad curated and distributed by "
                f"CZ CELLxGENE Discover in Collection: https://domain/collections/{collection.collection_id}"
            )
            adata = anndata.read_h5ad(f.name)
            self.assertEqual(adata.uns["citation"], citation_str)

    def test_process_validate_fail(self):
        """
        If the validation is not successful, the processing pipeline should:
        1. Set the processing status to INVALID
        2. Set a validation message accordingly
        """
        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        manifest = to_manifest(dropbox_uri)
        collection = self.generate_unpublished_collection()
        _, _ = self.business_logic.ingest_dataset(collection.version_id, dropbox_uri, None, None)

        # Set a mock failure for the schema validator
        self.schema_validator.validate_and_save_labels = Mock(
            return_value=(False, ["Validation error 1", "Validation error 2"], True)
        )

        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )

        pm = ProcessMain(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)

        for step_name in ["validate"]:
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
        self.assertEqual(status.validation_message, "Validation error 1\nValidation error 2")
