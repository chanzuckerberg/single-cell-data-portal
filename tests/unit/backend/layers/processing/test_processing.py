from unittest.mock import Mock, patch

from backend.layers.processing.downloader import Downloader
from backend.layers.processing.process import ProcessMain
from backend.layers.processing.process_cxg import ProcessCxg
from backend.layers.processing.process_download_validate import ProcessDownloadValidate
from backend.layers.processing.process_seurat import ProcessSeurat
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProviderInterface
from backend.layers.thirdparty.uri_provider import FileInfo, UriProvider
from backend.layers.common.entities import (
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetUploadStatus,
    DatasetValidationStatus,
)
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

    def upload_directory(self, src_dir: str, s3_uri: str):
        self.mock_s3_fs.append(s3_uri)

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
    def setUp(self):
        super().setUp()
        self.uri_provider = UriProvider()
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "local.h5ad"))
        self.s3_provider = MockS3Provider()
        self.schema_validator = SchemaValidatorProviderInterface()
        self.schema_validator.validate_and_save_labels = Mock(return_value=(True, [], True))
        self.downloader = Downloader(self.business_logic)
        self.downloader.download_file = Mock()

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
        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"

        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(collection.version_id, dropbox_uri, None)
        # This is where we're at when we start the SFN

        status = self.business_logic.get_dataset_status(dataset_version_id)
        # self.assertEqual(status.validation_status, DatasetValidationStatus.NA)
        self.assertIsNone(status.validation_status)
        self.assertEqual(status.processing_status, DatasetProcessingStatus.PENDING)
        self.assertEqual(status.upload_status, DatasetUploadStatus.WAITING)

        # TODO: ideally use a real h5ad so that
        with patch("backend.layers.processing.process_download_validate.ProcessDownloadValidate.extract_metadata"):
            pdv = ProcessDownloadValidate(
                self.business_logic, self.uri_provider, self.s3_provider, self.downloader, self.schema_validator
            )
            pdv.process(dataset_version_id, dropbox_uri, "fake_bucket_name")

            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.validation_status, DatasetValidationStatus.VALID)
            self.assertEqual(status.upload_status, DatasetUploadStatus.UPLOADED)
            self.assertEqual(status.h5ad_status, DatasetConversionStatus.UPLOADED)

            # Verify that both the original (raw.h5ad) and the labeled (local.h5ad) files are there
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/raw.h5ad"))
            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"))

            artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
            self.assertEqual(2, len(artifacts))
            artifact = artifacts[0]
            artifact.type = "H5AD"
            artifact.uri = f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"

    def test_process_seurat_success(self):
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(collection.version_id, "nothing", None)

        with patch("backend.layers.processing.process_seurat.ProcessSeurat.make_seurat") as mock:
            mock.return_value = "local.rds"
            ps = ProcessSeurat(self.business_logic, self.uri_provider, self.s3_provider)
            ps.process(dataset_version_id, "fake_bucket_name")

            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.rds_status, DatasetConversionStatus.UPLOADED)

            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.rds"))

            artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
            self.assertEqual(1, len(artifacts))
            artifact = artifacts[0]
            artifact.type = "RDS"
            artifact.uri = f"s3://fake_bucket_name/{dataset_version_id.id}/local.rds"

    def test_process_cxg_success(self):
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(collection.version_id, "nothing", None)

        with patch("backend.layers.processing.process_cxg.ProcessCxg.make_cxg") as mock:
            mock.return_value = "local.cxg"
            ps = ProcessCxg(self.business_logic, self.uri_provider, self.s3_provider)
            ps.process(dataset_version_id, "fake_bucket_name", "fake_cxg_bucket")

            status = self.business_logic.get_dataset_status(dataset_version_id)
            self.assertEqual(status.cxg_status, DatasetConversionStatus.UPLOADED)

            self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/"))

            artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
            self.assertEqual(1, len(artifacts))
            artifact = artifacts[0]
            artifact.type = "CXG"
            artifact.uri = f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/"

    @patch("backend.layers.processing.process_download_validate.ProcessDownloadValidate.extract_metadata")
    @patch("backend.layers.processing.process_seurat.ProcessSeurat.make_seurat")
    @patch("backend.layers.processing.process_cxg.ProcessCxg.make_cxg")
    def test_process_all(self, mock_cxg, mock_seurat, mock_h5ad):
        mock_seurat.return_value = "local.rds"
        mock_cxg.return_value = "local.cxg"
        dropbox_uri = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        collection = self.generate_unpublished_collection()
        dataset_version_id, dataset_id = self.business_logic.ingest_dataset(collection.version_id, dropbox_uri, None)

        pm = ProcessMain(
            self.business_logic, self.uri_provider, self.s3_provider, self.downloader, self.schema_validator
        )
        for step_name in ["download-validate", "cxg", "seurat"]:
            pm.process(dataset_version_id, step_name, dropbox_uri, "fake_bucket_name", "fake_cxg_bucket")

        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/raw.h5ad"))
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.h5ad"))
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_bucket_name/{dataset_version_id.id}/local.rds"))
        self.assertTrue(self.s3_provider.uri_exists(f"s3://fake_cxg_bucket/{dataset_version_id.id}.cxg/"))

        status = self.business_logic.get_dataset_status(dataset_version_id)
        self.assertEqual(status.cxg_status, DatasetConversionStatus.UPLOADED)
        self.assertEqual(status.rds_status, DatasetConversionStatus.UPLOADED)
        self.assertEqual(status.h5ad_status, DatasetConversionStatus.UPLOADED)
        self.assertEqual(status.validation_status, DatasetValidationStatus.VALID)
        self.assertEqual(status.upload_status, DatasetUploadStatus.UPLOADED)
        self.assertEqual(status.processing_status, DatasetProcessingStatus.PENDING)
        # TODO: DatasetProcessingStatus.SUCCESS is set by a lambda that also needs to be modified. It should belong here

        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_version_id))
        self.assertEqual(4, len(artifacts))

    def test_process_download_validate_fail(self):
        pass
