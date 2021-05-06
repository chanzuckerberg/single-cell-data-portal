from backend.corpora.common.corpora_orm import (
    DatasetArtifactFileType,
    DatasetArtifactType,
    DbDatasetProcessingStatus,
    UploadStatus,
)
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import processing_status_updater
from backend.corpora.lambdas.upload_failures.upload import update_dataset_processing_status_to_failed
from tests.unit.backend.corpora.common.entities.datasets import TestDataset
from tests.unit.backend.utils import BogusProcessingStatusParams, BogusDatasetParams


class TestUpdateDataset(TestDataset):
    def test__update__ok(self):
        artifact_params = dict(
            filename="filename_1",
            filetype=DatasetArtifactFileType.H5AD,
            type=DatasetArtifactType.ORIGINAL,
            user_submitted=True,
            s3_uri="some_uri",
        )
        deployment_directory_params = dict(url="test_url")
        processing_status = BogusProcessingStatusParams.get()
        dataset_params = BogusDatasetParams.get()

        dataset = Dataset.create(
            self.session,
            **dataset_params,
            artifacts=[artifact_params],
            deployment_directories=[deployment_directory_params],
            processing_status=processing_status,
        )

        new_artifact_params = dict(
            filename="a_different_filename",
            filetype=DatasetArtifactFileType.LOOM,
            type=DatasetArtifactType.ORIGINAL,
            user_submitted=False,
            s3_uri="a_different_uri",
        )
        new_processing_status = BogusProcessingStatusParams.get(upload_progress=7 / 9)

        dataset.update(
            artifacts=[new_artifact_params],
            deployment_directories=[deployment_directory_params],
            processing_status=new_processing_status,
            sex=["other"],
        )
        self.session.expire_all()
        actual_dataset = Dataset.get(self.session, dataset.id)
        self.assertEqual(actual_dataset.artifacts[0].filename, "a_different_filename")
        self.assertEqual(actual_dataset.sex, ["other"])
        self.assertEqual(actual_dataset.processing_status.upload_progress, 7 / 9)

    def test__update_processing_status__ok(self):
        dataset = Dataset.get(self.session, self.uuid)
        status = {
            DbDatasetProcessingStatus.upload_progress: 0,
            DbDatasetProcessingStatus.upload_status: UploadStatus.WAITING,
        }

        processing_status_updater(self.session, dataset.processing_status.id, status)

        dataset = Dataset.get(self.session, self.uuid)
        self.assertEqual(dataset.processing_status.upload_status, UploadStatus.WAITING)
        update_dataset_processing_status_to_failed(self.uuid)
        self.session.expire_all()

        dataset = Dataset.get(self.session, self.uuid)
        self.assertEqual(dataset.processing_status.upload_status, UploadStatus.FAILED)

    def test__update_processing_status__no_dataset__ok(self):
        update_dataset_processing_status_to_failed("fake_uuid")
