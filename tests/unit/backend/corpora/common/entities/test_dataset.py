import typing
import unittest

from backend.corpora.common.corpora_orm import (
    DbDatasetArtifact,
    DbDeploymentDirectory,
    DatasetArtifactType,
    DatasetArtifactFileType,
    CollectionVisibility,
    DbDataset,
    DbCollection,
    Base,
)
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.utils.db_utils import DbUtils
from tests.unit.backend.utils import BogusDatasetParams


class TestDataset(unittest.TestCase):
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

    def test__get__does_not_exist(self):
        non_existent_key = "non_existent_id"
        self.assertEqual(Dataset.get(non_existent_key), None)

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

        deployment_directory_params = dict(url="test_url")

        dataset_params = BogusDatasetParams.get()

        for i in range(3):
            with self.subTest(i):
                dataset = Dataset.create(
                    **dataset_params,
                    artifacts=[artifact_params] * i,
                    deployment_directories=[deployment_directory_params] * i,
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

    def test__list__ok(self):
        generate = 2
        generated_ids = [Dataset.create(**BogusDatasetParams.get()).id for _ in range(generate)]
        dataset = Dataset.list()
        self.assertTrue(set(generated_ids).issubset([d.id for d in dataset]))

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

    def test__get_asset__ok(self):
        dataset = Dataset.get(self.uuid)
        expected_asset_id = "test_dataset_artifact_id"
        asset = dataset.get_asset("test_dataset_artifact_id")
        self.assertEqual(expected_asset_id, asset.id)

    def test__get_asset__not_found(self):
        dataset = Dataset.get(self.uuid)
        asset = dataset.get_asset("fake_asset")
        self.assertIsNone(asset)

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
