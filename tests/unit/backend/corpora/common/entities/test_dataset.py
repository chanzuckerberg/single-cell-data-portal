import unittest

from backend.corpora.common.corpora_orm import (
    DbContributor,
    DbDatasetArtifact,
    DbDeploymentDirectory,
    DatasetArtifactType,
    DatasetArtifactFileType,
)
from backend.corpora.common.entities.dataset import Dataset
from tests.unit.backend.corpora.common.entities.utils import get_ids


class DatasetParams:
    @classmethod
    def get(cls):
        return dict(
            name=f"create_dataset",
            organism="organism",
            organism_ontology="123",
            tissue="tissue",
            tissue_ontology="123",
            assay="assay",
            assay_ontology="123",
            disease="diseas",
            disease_ontology="123",
            sex="F",
            ethnicity="ethnicity",
            ethnicity_ontology="123",
            source_data_location="location",
            preprint_doi="preprint",
            publication_doi="publication",
        )


class TestDataset(unittest.TestCase):
    def setUp(self):
        self.uuid = "test_dataset_id"

    def test__get__ok(self):
        dataset = Dataset.get(self.uuid)
        self.assertEqual(dataset.id, self.uuid)
        self.assertEqual(dataset.assay, "test_assay")

        # Verify Artifact relationship
        self.assertIsInstance(dataset.artifacts[0], DbDatasetArtifact)
        self.assertEqual(dataset.artifacts[0].id, "test_dataset_artifact_id")

        # Verify Deployment Directory relationship
        self.assertIsInstance(dataset.deployment_directories[0], DbDeploymentDirectory)
        self.assertEqual(dataset.deployment_directories[0].id, "test_deployment_directory_id")

        # Verify Contributor Relationship
        self.assertIsInstance(dataset.contributors[0], DbContributor)
        self.assertEqual(dataset.contributors[0].id, "test_contributor_id")

    def test__get__does_not_exist(self):
        non_existent_key = "non_existent_id"
        self.assertEqual(Dataset.get(non_existent_key), None)

    def test__create__ok(self):
        """
        Create a dataset with a variable number of artifacts, contributors, and deployment_directories
        """
        artifact_params = dict(
            filename="filename_1",
            filetype=DatasetArtifactFileType.H5AD,
            type=DatasetArtifactType.ORIGINAL,
            user_submitted=True,
            s3_uri="some_uri",
        )

        deployment_directory_params = dict(environment="test", url="test_url")

        contributor_params = dict(name="bob", institution="school", email="some@email.com")

        dataset_params = DatasetParams.get()

        for i in range(3):
            with self.subTest(i):
                dataset = Dataset.create(
                    **dataset_params,
                    artifacts=[artifact_params] * i,
                    contributors=[contributor_params] * i,
                    deployment_directories=[deployment_directory_params] * i,
                )

                dataset_id = dataset.id
                artifact_ids = get_ids(dataset.artifacts)
                deployment_directory_ids = get_ids(dataset.deployment_directories)
                contributor_ids = get_ids(dataset.contributors)

                # Expire all local objects and retrieve them from the DB to make sure the transactions went through.
                Dataset.db.session.expire_all()

                dataset = Dataset.get(dataset_id)
                self.assertIsNotNone(dataset)
                self.assertEqual(dataset_id, dataset.id)
                self.assertEqual(artifact_ids, get_ids(dataset.artifacts))
                self.assertEqual(deployment_directory_ids, get_ids(dataset.deployment_directories))
                self.assertEqual(contributor_ids, get_ids(dataset.contributors))

    def test__create_ids__ok(self):
        """
        Creating a dataset with ids in connect attributes. A new id is generated even if id is provided.
        """
        dataset_params = DatasetParams.get()
        dataset = Dataset.create(
            **dataset_params,
            artifacts=[{"id": "test_dataset_artifact_id"}],
            contributors=[{"id": "test_contributor_id"}],
            deployment_directories=[{"id": "test_deployment_directory_id"}],
        )

        dataset_id = dataset.id
        artifact_ids = get_ids(dataset.artifacts)
        deployment_directory_ids = get_ids(dataset.deployment_directories)
        contributor_ids = get_ids(dataset.contributors)

        # Expire all local objects and retrieve them from the DB to make sure the transactions went through.
        Dataset.db.session.expire_all()

        dataset = Dataset.get(dataset_id)
        self.assertIsNotNone(dataset)
        self.assertEqual(dataset_id, dataset.id)

        self.assertEqual(artifact_ids, get_ids(dataset.artifacts))
        self.assertNotEqual(["test_dataset_artifact_id"], get_ids(dataset.artifacts))

        self.assertEqual(deployment_directory_ids, get_ids(dataset.deployment_directories))
        self.assertNotEqual(["test_dataset_artifact_id"], get_ids(dataset.deployment_directories))

        self.assertEqual(contributor_ids, get_ids(dataset.contributors))
        self.assertEqual(["test_contributor_id"], get_ids(dataset.contributors))

    def test__list__ok(self):
        generate = 5

        for _ in range(generate):
            Dataset.create(**DatasetParams.get())

        datasets = Dataset.list()
        self.assertGreaterEqual(len(datasets), generate)
        self.assertTrue(all([isinstance(i, Dataset) for i in datasets]))
