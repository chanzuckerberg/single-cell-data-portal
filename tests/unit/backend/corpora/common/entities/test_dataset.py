import unittest

from backend.corpora.common.corpora_orm import (
    DbContributor,
    DbDatasetArtifact,
    DbDeploymentDirectory,
    DatasetArtifactType,
    DatasetArtifactFileType,
)
from backend.corpora.common.entities.dataset import Dataset


class DatasetParams:
    @classmethod
    def get(cls):
        return dict(
            name="create_dataset",
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

                expected_dataset_id = dataset.id
                expected_artifacts = dataset.artifacts
                expected_deployment_directories = dataset.deployment_directories
                expected_contributors = dataset.contributors

                # Expire all local objects and retrieve them from the DB to make sure the transactions went through.
                Dataset.db.session.expire_all()

                actual_dataset = Dataset.get(expected_dataset_id)
                self.assertEqual(expected_dataset_id, actual_dataset.id)
                self.assertCountEqual(expected_artifacts, actual_dataset.artifacts)
                self.assertCountEqual(expected_deployment_directories, actual_dataset.deployment_directories)
                self.assertCountEqual(expected_contributors, actual_dataset.contributors)
