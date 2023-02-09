import logging
from unittest import mock
from unittest.mock import patch

from backend.common.corpora_orm import ConversionStatus, DatasetArtifactFileType
from backend.portal.pipeline.processing.exceptions import ConversionFailed
from backend.portal.pipeline.processing.process_seurat import process
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestProcessSeurat(CorporaTestCaseUsingMockAWS):
    @patch("backend.portal.pipeline.processing.process_seurat.download_from_s3")
    @patch("backend.portal.pipeline.processing.process_seurat.make_seurat")
    def test_process_with_seurat_conversion_failures(self, mock_seurat, mock_download_from_s3):
        mock_seurat.side_effect = RuntimeError("seurat conversion failed")
        # given
        dataset = self.generate_dataset(self.session)
        self.generate_artifact(self.session, dataset.id, DatasetArtifactFileType.H5AD)
        # then
        with self.assertRaises(ConversionFailed):
            process(dataset.id, self.bucket.name)

    @patch("backend.portal.pipeline.processing.process_seurat.make_seurat")
    def test__process_skips_seurat_conversion_when_unconvertible_dataset_detected(self, mock_make_seurat):
        # given
        dataset = self.generate_dataset(self.session, processing_status={"rds_status": ConversionStatus.SKIPPED})

        # when
        process(dataset.id, mock.ANY)

        # then
        mock_make_seurat.assert_not_called()

    @patch("backend.portal.pipeline.processing.process_seurat.make_seurat")
    @patch("backend.portal.pipeline.processing.process_seurat.download_from_s3")
    @patch("backend.portal.pipeline.processing.process_seurat.create_artifact")
    def test__process_runs_seurat_conversion_when_convertible_dataset_detected_1(
        self,
        mock_download_from_s3,
        mock_create_artifact,
        mock_make_seurat,
    ):
        # given
        dataset = self.generate_dataset(self.session)
        self.generate_artifact(self.session, dataset.id, DatasetArtifactFileType.H5AD)
        # when
        process(dataset.id, self.bucket.name)

        # then
        mock_download_from_s3.assert_called()
        mock_make_seurat.assert_called()
        mock_create_artifact.assert_called()

    @patch("backend.portal.pipeline.processing.process_seurat.make_seurat")
    @patch("backend.portal.pipeline.processing.process_seurat.download_from_s3")
    @patch("backend.portal.pipeline.processing.process_seurat.replace_artifact")
    def test__process_runs_seurat_conversion_when_convertible_dataset_detected_2(
        self,
        mock_replace_artifact,
        mock_download_from_s3,
        mock_make_seurat,
    ):
        # given
        dataset = self.generate_dataset(self.session)
        self.generate_artifact(self.session, dataset.id, DatasetArtifactFileType.H5AD)
        self.generate_artifact(self.session, dataset.id, DatasetArtifactFileType.RDS)
        # when
        with self.assertLogs("dataset_processing", logging.WARNING) as logs:
            process(dataset.id, self.bucket.name)

        # then
        self.assertIn("will replace the S3 file only", logs.output[0])
        mock_download_from_s3.assert_called()
        mock_make_seurat.assert_called()
        mock_replace_artifact.assert_called()

    @patch("backend.portal.pipeline.processing.process_seurat.make_seurat")
    @patch("backend.portal.pipeline.processing.process_seurat.download_from_s3")
    @patch("backend.portal.pipeline.processing.process_seurat.create_artifact")
    def test__process_runs_seurat_conversion_when_convertible_dataset_detected_3(
        self,
        mock_download_from_s3,
        mock_create_artifact,
        mock_make_seurat,
    ):
        # given
        dataset = self.generate_dataset(self.session)
        artifact = self.generate_artifact(
            self.session,
            dataset.id,
            DatasetArtifactFileType.H5AD,
        )
        artifact.update(s3_uri="bucket/path/uuid/file")

        # when
        with self.assertLogs("dataset_processing", logging.WARNING) as logs:
            process(dataset.id, self.bucket.name)

        # then
        self.assertIn("creating a new artifact", logs.output[0])
        mock_download_from_s3.assert_called()
        mock_make_seurat.assert_called()
        mock_create_artifact.assert_called()

    @patch("backend.portal.pipeline.processing.process_seurat.make_seurat")
    @patch("backend.portal.pipeline.processing.process_seurat.download_from_s3")
    @patch("backend.portal.pipeline.processing.process_seurat.replace_artifact")
    def test__process_runs_seurat_conversion_when_convertible_dataset_detected_4(
        self,
        mock_replace_artifact,
        mock_download_from_s3,
        mock_make_seurat,
    ):
        # given
        dataset = self.generate_dataset(self.session)
        artifact = self.generate_artifact(
            self.session,
            dataset.id,
            DatasetArtifactFileType.H5AD,
        )
        artifact.update(s3_uri="bucket/path/uuid/file")
        self.generate_artifact(self.session, dataset.id, DatasetArtifactFileType.RDS)

        # when
        with self.assertLogs("dataset_processing", logging.WARNING) as logs:
            process(dataset.id, self.bucket.name)

        # then
        self.assertIn("replacing it", logs.output[0])
        mock_download_from_s3.assert_called()
        mock_make_seurat.assert_called()
        mock_replace_artifact.assert_called()
