from backend.corpora.common.corpora_orm import (
    DatasetArtifactFileType,
    DatasetArtifactType,
    UploadStatus,
    ValidationStatus,
)
from backend.corpora.common.entities import Dataset
from tests.unit.backend.corpora.common.entities.datasets import TestDataset
from tests.unit.backend.utils import BogusDatasetParams


class TestCreateDataset(TestDataset):
    def test__create__ok(self):

        """
        Create a dataset with a variable number of artifacts, and deployment_directories
        """
        artifact_params = dict(
            filename="filename_1",
            filetype=DatasetArtifactFileType.H5AD,
            type=DatasetArtifactType.ORIGINAL,
            user_submitted=True,
            s3_uri="some_uri",
        )

        for i in range(3):
            with self.subTest(i):
                dataset = self.create_dataset_with_artifacts(
                    artifact_count=i, deployment_dir_count=i, artifact_params=artifact_params
                )

                expected_dataset_id = dataset.id
                expected_artifacts = [art.to_dict() for art in dataset.artifacts]
                expected_deployment_directories = [dep.to_dict() for dep in dataset.deployment_directories]

                # Expire all local objects and retrieve them from the DB to make sure the transactions went through.
                self.session.expire_all()

                actual_dataset = Dataset.get(self.session, expected_dataset_id)
                actual_artifacts = [art.to_dict() for art in actual_dataset.artifacts]
                actual_deployment_directories = [dep.to_dict() for dep in actual_dataset.deployment_directories]
                self.assertEqual(expected_dataset_id, actual_dataset.id)
                self.assertCountEqual(expected_artifacts, actual_artifacts)
                self.assertCountEqual(expected_deployment_directories, actual_deployment_directories)

    def test__create_with_processing_status(self):
        """
        Create a dataset with a processing status
        """

        dataset_params = BogusDatasetParams.get()

        dataset = Dataset.create(
            self.session,
            **dataset_params,
            processing_status={
                "upload_progress": 9 / 13,
                "upload_status": UploadStatus.UPLOADING,
                "validation_status": ValidationStatus.NA,
            },
        )
        expected_dataset_id = dataset.id

        # Expire all local objects and retrieve them from the DB to make sure the transactions went through.
        self.session.expire_all()

        actual_dataset = Dataset.get(self.session, expected_dataset_id)
        self.assertAlmostEqual(actual_dataset.processing_status.upload_progress, 9 / 13)
        self.assertEqual(actual_dataset.processing_status.upload_status, UploadStatus.UPLOADING)
        self.assertEqual(actual_dataset.processing_status.validation_status, ValidationStatus.NA)
