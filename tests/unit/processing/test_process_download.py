from unittest.mock import Mock

import pytest

from backend.layers.common.entities import DatasetArtifactType, DatasetProcessingStatus, DatasetUploadStatus
from backend.layers.processing.process_download import ProcessDownload
from tests.unit.processing.base_processing_test import BaseProcessingTest


class TestProcessDownload(BaseProcessingTest):
    def test_process_download_success(self):
        """
        ProcessValidate should:
        1. Download the h5ad artifact
        2. Set upload status to UPLOADED
        3. upload the original file to S3
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

        pdv = ProcessDownload(self.business_logic, self.uri_provider, self.s3_provider)
        pdv.process(dataset_version_id, dropbox_uri, "fake_bucket_name", "fake_sfn_task_token")

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.upload_status, DatasetUploadStatus.UPLOADED)

        # Verify that both the original (raw.h5ad) and the labeled (local.h5ad) files are there
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/raw.h5ad"))

        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
        self.assertEqual(1, len(artifacts))
        artifact = artifacts[0]
        artifact.type = DatasetArtifactType.RAW_H5AD
        artifact.uri = f"s3://fake_bucket_name/{dataset_version_id.id}/raw.h5ad"

    def test_download_from_s3_uri(self):
        """
        Call process download using an s3 uri
        """

        s3_uri = "s3://fake_bucket_name/fake_key/fake_file.h5ad"
        pdv = ProcessDownload(self.business_logic, self.uri_provider, self.s3_provider)
        pdv.download_from_s3 = Mock()

        assert pdv.download_from_source_uri(s3_uri, "fake_local_path") == "fake_local_path"
        pdv.download_from_s3.assert_called_once_with(
            bucket_name="fake_bucket_name", object_key="fake_key/fake_file.h5ad", local_filename="fake_local_path"
        )

    def test_download_unknown_uri(self):
        """
        Call process download using unknown
        """

        uri = "fake://fake_bucket_name/fake_key/fake_file.h5ad"
        pdv = ProcessDownload(self.business_logic, self.uri_provider, self.s3_provider)
        pdv.download_from_s3 = Mock()
        with pytest.raises(ValueError, match=f"Malformed source URI: {uri}"):
            pdv.download_from_source_uri(uri, "fake_local_path")
