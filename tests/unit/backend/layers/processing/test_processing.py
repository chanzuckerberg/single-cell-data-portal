import tempfile
import unittest
from unittest.mock import Mock, patch

from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProviderInterface
from backend.layers.processing.downloader import Downloader
from backend.layers.processing.process_download_validate import ProcessDownloadValidate
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProvider
from backend.layers.thirdparty.uri_provider import FileInfo, UriProvider
from backend.layers.common.entities import DatasetConversionStatus, DatasetProcessingStatus, DatasetUploadStatus, DatasetValidationStatus
from tests.unit.backend.layers.common.base_api_test import NewBaseTest


class MockS3Provider(S3Provider):
    """
    Simple S3 mock that mostly checks if paths are correct
    """

    def __init__(self) -> None:
        self.mock_s3_fs = []

    def upload_file(self, src_file: str, bucket_name: str, dst_file: str, extra_args: dict):
        url = f"s3://{bucket_name}/{dst_file}"
        self.mock_s3_fs.append(url)

    def download_file(self, bucket_name: str, object_key: str, local_filename: str):
        pass

    def file_exists(self, bucket_name: str, object_key: str):
        url = f"s3://{bucket_name}/{object_key}"
        return self.uri_exists(url)

    def uri_exists(self, uri: str):
        return uri in self.mock_s3_fs

    def is_empty(self):
        return len(self.mock_s3_fs) == 0


class ProcessingTest(NewBaseTest):


    def test_process_download_validate_success(self):
        """
        ProcessDownloadValidate should do things
        """
        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"

        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(collection.version_id, dropbox_uri, None)
        # This is where we're at when we start the SFN

        self.uri_provider = UriProvider()
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "local.h5ad"))
        self.s3_provider = MockS3Provider()
        schema_validator = SchemaValidatorProvider()
        schema_validator.validate_and_save_labels = Mock(return_value=(True, [], True))
        downloader = Downloader(self.business_logic)
        downloader.download_file = Mock()

        status = self.business_logic.get_dataset_status(dataset_version_id)
        # self.assertEqual(status.validation_status, DatasetValidationStatus.NA)
        self.assertIsNone(status.validation_status)
        self.assertEqual(status.processing_status, DatasetProcessingStatus.PENDING)
        self.assertEqual(status.upload_status, DatasetUploadStatus.WAITING)

        # TODO: ideally use a real h5ad so that 
        with patch("backend.layers.processing.process_download_validate.ProcessDownloadValidate.extract_metadata") as mock:
            pdv = ProcessDownloadValidate(self.business_logic, self.uri_provider, self.s3_provider, downloader, schema_validator)
            pdv.process(dataset_version_id, dropbox_uri, "fake_bucket_name")

            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.validation_status, DatasetValidationStatus.VALID)
            self.assertEqual(status.upload_status, DatasetUploadStatus.UPLOADED)
            self.assertEqual(status.h5ad_status, DatasetConversionStatus.UPLOADED)

            # Verify that both the original (raw.h5ad) and the labeled (local.h5ad) files are there
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/raw.h5ad"))
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"))



    def test_process_download_validate_fail(self):
        pass