import typing

from backend.corpora.common.corpora_orm import (
    DbDatasetArtifact,
    DbDatasetProcessingStatus,
    DbDeploymentDirectory,
    DatasetArtifactType,
    DatasetArtifactFileType,
    CollectionVisibility,
    UploadStatus,
    ValidationStatus,
    DbDataset,
    DbCollection,
    Base,
)
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.utils.db_utils import DbUtils, processing_status_updater
from backend.corpora.lambdas.upload_failures.upload import update_dataset_processing_status_to_failed
from tests.unit.backend.fixtures.generate_data_mixin import GenerateDataMixin
from tests.unit.backend.utils import BogusDatasetParams, BogusProcessingStatusParams
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestDataset(DataPortalTestCase, GenerateDataMixin):
    def setUp(self):
        self.uuid = "test_dataset_id"

    def test__get__ok(self):
        dataset = Dataset.get(self.uuid)
        self.assertEqual(dataset.id, self.uuid)
        self.assertEqual(len(dataset.assay), 1)
        self.assertDictEqual(dataset.assay[0], {"ontology_term_id": "test_obo", "label": "test_assay"})

        # Verify Artifact relationship
        self.assertIsInstance(dataset.artifacts[0], DbDatasetArtifact)
        self.assertEqual(dataset.artifacts[0].id, "test_dataset_artifact_id")

        # Verify Deployment Directory relationship
        self.assertIsInstance(dataset.deployment_directories[0], DbDeploymentDirectory)
        self.assertEqual(dataset.deployment_directories[0].id, "test_deployment_directory_id")

        # Verify Processing Status relationship
        self.assertIsInstance(dataset.processing_status, DbDatasetProcessingStatus)
        self.assertEqual(dataset.processing_status.id, "test_dataset_processing_status_id")
        self.assertEqual(dataset.processing_status.dataset_id, "test_dataset_id")

    def test__get__does_not_exist(self):
        non_existent_key = "non_existent_id"
        self.assertEqual(Dataset.get(non_existent_key), None)

    def create_dataset_with_artifacts(self, artifact_count=1, deployment_dir_count=1, artifact_params=None):
        """
        Create a dataset with a variable number of artifacts, and deployment_directories
        """
        if not artifact_params:
            artifact_params = dict(
                filename="filename_x",
                filetype=DatasetArtifactFileType.H5AD,
                type=DatasetArtifactType.ORIGINAL,
                user_submitted=True,
                s3_uri="some_uri",
            )

        deployment_directory_params = dict(url="test_url")

        dataset_params = BogusDatasetParams.get()
        dataset = self.generate_dataset(
            **dataset_params,
            artifacts=[artifact_params] * artifact_count,
            deployment_directories=[deployment_directory_params] * deployment_dir_count,
            processing_status={
                "upload_progress": 9 / 13,
                "upload_status": UploadStatus.UPLOADING,
                "validation_status": ValidationStatus.NA,
            },
        )
        return dataset

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
                Dataset.db.session.expire_all()

                actual_dataset = Dataset.get(expected_dataset_id)
                actual_artifacts = [art.to_dict() for art in actual_dataset.artifacts]
                actual_deployment_directories = [dep.to_dict() for dep in actual_dataset.deployment_directories]
                self.assertEqual(expected_dataset_id, actual_dataset.id)
                self.assertCountEqual(expected_artifacts, actual_artifacts)
                self.assertCountEqual(expected_deployment_directories, actual_deployment_directories)

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
        Dataset.db.session.expire_all()
        actual_dataset = Dataset.get(dataset.id)
        self.assertEqual(actual_dataset.artifacts[0].filename, "a_different_filename")
        self.assertEqual(actual_dataset.sex, ["other"])
        self.assertEqual(actual_dataset.processing_status.upload_progress, 7 / 9)

    def test__list__ok(self):
        generate = 2
        generated_ids = [Dataset.create(**BogusDatasetParams.get()).id for _ in range(generate)]
        dataset = Dataset.list()
        self.assertTrue(set(generated_ids).issubset([d.id for d in dataset]))

    def test__create_with_processing_status(self):
        """
        Create a dataset with a processing status
        """

        dataset_params = BogusDatasetParams.get()

        dataset = Dataset.create(
            **dataset_params,
            processing_status={
                "upload_progress": 9 / 13,
                "upload_status": UploadStatus.UPLOADING,
                "validation_status": ValidationStatus.NA,
            },
        )
        expected_dataset_id = dataset.id

        # Expire all local objects and retrieve them from the DB to make sure the transactions went through.
        Dataset.db.session.expire_all()

        actual_dataset = Dataset.get(expected_dataset_id)
        self.assertAlmostEqual(actual_dataset.processing_status.upload_progress, 9 / 13)
        self.assertEqual(actual_dataset.processing_status.upload_status, UploadStatus.UPLOADING)
        self.assertEqual(actual_dataset.processing_status.validation_status, ValidationStatus.NA)

    def test__cascade_delete_dataset__ok(self):
        # Create the dataset
        test_dataset = Dataset.create(
            **BogusDatasetParams.get(
                collection_id="test_collection_id",
                collection_visibility=CollectionVisibility.PUBLIC.name,
                artifacts=[{}],
                deployment_directories=[{}],
            )
        )
        test_dataset_ids = [(test_dataset.id, DbDataset)]
        test_artifact_ids = [(art.id, DbDatasetArtifact) for art in test_dataset.artifacts]
        test_deployed_directory_ids = [(dep.id, DbDeploymentDirectory) for dep in test_dataset.deployment_directories]
        test_collection_ids = [(("test_collection_id", CollectionVisibility.PUBLIC.name), DbCollection)]

        with self.subTest("verify everything exists"):
            expected_exists = test_collection_ids + test_dataset_ids + test_artifact_ids + test_deployed_directory_ids
            self.assertRowsExist(expected_exists)

        # Delete the dataset
        test_dataset.delete()

        with self.subTest("Verify Deletion"):
            expected_deleted = test_dataset_ids + test_artifact_ids + test_deployed_directory_ids
            expected_exists = test_collection_ids
            self.assertRowsDeleted(expected_deleted)
            self.assertRowsExist(expected_exists)

    def test__update_processing_status__ok(self):
        dataset = Dataset.get(self.uuid)
        status = {
            DbDatasetProcessingStatus.upload_progress: 0,
            DbDatasetProcessingStatus.upload_status: UploadStatus.WAITING,
        }

        processing_status_updater(dataset.processing_status.id, status)

        dataset = Dataset.get(self.uuid)
        self.assertEqual(dataset.processing_status.upload_status, UploadStatus.WAITING)
        update_dataset_processing_status_to_failed(self.uuid)
        Dataset.db.session.expire_all()

        dataset = Dataset.get(self.uuid)
        self.assertEqual(dataset.processing_status.upload_status, UploadStatus.FAILED)

    def test__update_processing_status__no_dataset__ok(self):
        update_dataset_processing_status_to_failed("fake_uuid")

    def test__get_asset__ok(self):
        dataset = Dataset.get(self.uuid)
        expected_asset_id = "test_dataset_artifact_id"
        asset = dataset.get_asset("test_dataset_artifact_id")
        self.assertEqual(expected_asset_id, asset.id)

    def test__get_asset__not_found(self):
        dataset = Dataset.get(self.uuid)
        asset = dataset.get_asset("fake_asset")
        self.assertIsNone(asset)

    def test__tombstone_dataset_and_delete_child_objects(self):
        dataset = self.create_dataset_with_artifacts(artifact_count=3, deployment_dir_count=2)
        self.assertEqual(dataset.processing_status.upload_status, UploadStatus.UPLOADING)
        self.assertEqual(len(dataset.artifacts), 3)
        self.assertEqual(len(dataset.deployment_directories), 2)
        self.assertFalse(dataset.tombstone)

        dataset.tombstone_dataset_and_delete_child_objects()
        self.assertEqual(len(dataset.artifacts), 0)
        self.assertEqual(len(dataset.deployment_directories), 0)
        self.assertTrue(dataset.tombstone)
        self.assertIsNone(dataset.processing_status)

    def assertRowsDeleted(self, tests: typing.List[typing.Tuple[str, Base]]):
        """
        Verify if rows have been deleted from the database.
        :param tests: a list of tuples with (primary_key, table)
        """
        db = DbUtils()
        db.session.expire_all()
        for p_key, table in tests:
            if len(p_key) == 2:
                # handle the special case for collections with a composite primary key
                actual = db.query([table], [table.id == p_key[0], table.visibility == p_key[1]])
            else:
                actual = db.query([table], [table.id == p_key])
            self.assertFalse(actual, f"Row not deleted {table.__name__}:{p_key}")

    def assertRowsExist(self, tests):
        """
        Verify if rows exist in the database.
        :param tests: a list of tuples with (primary_key, table)
        """
        db = DbUtils()
        db.session.expire_all()
        for p_key, table in tests:
            if len(p_key) == 2:
                # handle the special case for collections with a composite primary key
                actual = db.query([table], [table.id == p_key[0], table.visibility == p_key[1]])
            else:
                actual = db.query([table], [table.id == p_key])
            self.assertTrue(actual, f"Row does not exist {table.__name__}:{p_key}")
