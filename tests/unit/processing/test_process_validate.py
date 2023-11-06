from unittest.mock import MagicMock, Mock, patch

from backend.layers.common.entities import (
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersionId,
    Link,
)
from backend.layers.processing.process import ProcessMain
from backend.layers.processing.process_validate import ProcessValidate
from tests.unit.processing.base_processing_test import BaseProcessingTest

mock_config_attr = {
    "collections_base_url": "http://collections",
    "dataset_assets_base_url": "http://domain",
    "upload_max_file_size_gb": 1,
    "schema_4_feature_flag": "True",
}


def mock_config_fn(name):
    return mock_config_attr[name]


class ProcessingTest(BaseProcessingTest):
    @patch("scanpy.read_h5ad")
    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    def test_process_download_validate_success(self, mock_config, mock_read_h5ad):
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

        collection = self.generate_unpublished_collection(
            links=[Link(name=None, type="DOI", uri="http://doi.org/12.2345")]
        )
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )
        # This is where we're at when we start the SFN

        status = self.business_logic.get_dataset_status(dataset_version_id)
        # self.assertEqual(status.validation_status, DatasetValidationStatus.NA)
        self.assertIsNone(status.validation_status)
        self.assertEqual(status.processing_status, DatasetProcessingStatus.INITIALIZED)
        self.assertEqual(status.upload_status, DatasetUploadStatus.WAITING)

        mock_read_h5ad.return_value = MagicMock(uns=dict())

        # TODO: ideally use a real h5ad so that
        with patch("backend.layers.processing.process_validate.ProcessValidate.extract_metadata"):
            pdv = ProcessValidate(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
            pdv.process(collection.version_id, dataset_version_id, "fake_bucket_name", "fake_datasets_bucket")
            citation_str = (
                f"Publication: http://doi.org/12.2345 "
                f"Dataset Version: http://domain/{dataset_version_id}.h5ad curated and distributed by "
                f"CZ CELLxGENE Discover in Collection: http://collections/{collection.version_id}"
            )
            self.assertEqual(mock_read_h5ad.return_value.uns["citation"], citation_str)
            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.validation_status, DatasetValidationStatus.VALID)
            self.assertEqual(status.h5ad_status, DatasetConversionStatus.UPLOADED)

            # Verify that both the original (raw.h5ad) and the labeled (local.h5ad) files are there
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"))
            # Verify that the labeled file is uploaded to the datasets bucket
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_datasets_bucket/{dataset_version_id.id}.h5ad"))

            artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
            self.assertEqual(1, len(artifacts))
            artifact = artifacts[0]
            artifact.type = DatasetArtifactType.H5AD
            artifact.uri = f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"

    @patch("scanpy.read_h5ad")
    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    def test_populate_dataset_citation__no_publication_doi(self, mock_config, mock_read_h5ad):
        mock_read_h5ad.return_value = MagicMock(uns=dict())
        collection = self.generate_unpublished_collection()

        pdv = ProcessValidate(self.business_logic, self.uri_provider, self.s3_provider, self.schema_validator)
        dataset_version_id = DatasetVersionId()
        pdv.populate_dataset_citation(collection.version_id, dataset_version_id, "")
        citation_str = (
            f"Dataset Version: http://domain/{dataset_version_id}.h5ad curated and distributed by "
            f"CZ CELLxGENE Discover in Collection: http://collections/{collection.version_id}"
        )
        self.assertEqual(mock_read_h5ad.return_value.uns["citation"], citation_str)

    def test_process_validate_fail(self):
        """
        If the validation is not successful, the processing pipeline should:
        1. Set the processing status to INVALID
        2. Set a validation message accordingly
        """
        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(
            collection.version_id, dropbox_uri, None, None
        )

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
                dropbox_uri,
                "fake_bucket_name",
                "fake_datasets_bucket",
                "fake_cxg_bucket",
            )

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.validation_status, DatasetValidationStatus.INVALID)
        self.assertEqual(status.validation_message, "Validation error 1\nValidation error 2")
