import tempfile

import filecmp
import os
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
from backend.corpora.common.entities.dataset import Dataset, get_cxg_bucket_path
from backend.corpora.common.entities.geneset import GenesetDatasetLink, Geneset
from backend.corpora.common.utils.db_session import processing_status_updater
from backend.corpora.lambdas.upload_failures.upload import update_dataset_processing_status_to_failed
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from tests.unit.backend.utils import BogusDatasetParams, BogusProcessingStatusParams


class TestDataset(CorporaTestCaseUsingMockAWS):
    def setUp(self):
        super().setUp()
        self.uuid = "test_dataset_id"
        self.bucket_name = self.CORPORA_TEST_CONFIG["bucket_name"]

    def test__get__ok(self):
        dataset = Dataset.get(self.session, self.uuid)
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
        self.assertEqual(Dataset.get(self.session, non_existent_key), None)

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
            self.session,
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
                self.session.expire_all()

                actual_dataset = Dataset.get(self.session, expected_dataset_id)
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

    def test__list__ok(self):
        generate = 2
        generated_ids = [Dataset.create(self.session, **BogusDatasetParams.get()).id for _ in range(generate)]
        dataset = Dataset.list(self.session)
        self.assertTrue(set(generated_ids).issubset([d.id for d in dataset]))

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

    def test__cascade_delete_dataset__ok(self):
        # Create the dataset
        test_dataset = Dataset.create(
            self.session,
            **BogusDatasetParams.get(
                collection_id="test_collection_id",
                collection_visibility=CollectionVisibility.PUBLIC.name,
                artifacts=[{}],
                deployment_directories=[{}],
            ),
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

    def test__get_asset__ok(self):
        dataset = Dataset.get(self.session, self.uuid)
        expected_asset_id = "test_dataset_artifact_id"
        asset = dataset.get_asset("test_dataset_artifact_id")
        self.assertEqual(expected_asset_id, asset.id)

    def test__get_asset__not_found(self):
        dataset = Dataset.get(self.session, self.uuid)
        asset = dataset.get_asset("fake_asset")
        self.assertIsNone(asset)

    def test__tombstone_dataset_and_delete_child_objects(self):
        dataset = self.create_dataset_with_artifacts(artifact_count=3, deployment_dir_count=2)
        geneset = self.generate_geneset(self.session)
        GenesetDatasetLink.create(self.session, geneset.id, dataset.id)
        self.assertEqual(dataset.processing_status.upload_status, UploadStatus.UPLOADING)
        self.assertEqual(len(dataset.artifacts), 3)
        self.assertEqual(len(dataset.deployment_directories), 2)
        self.assertEqual(len(dataset.genesets), 1)
        self.assertFalse(dataset.tombstone)

        dataset.tombstone_dataset_and_delete_child_objects()
        self.assertEqual(len(dataset.artifacts), 0)
        self.assertEqual(len(dataset.deployment_directories), 0)
        self.assertEqual(len(dataset.genesets), 0)
        self.assertTrue(dataset.tombstone)
        self.assertIsNone(dataset.processing_status)
        self.assertIsNotNone(Geneset.get(self.session, geneset.id))

    def test__deletes_assets_from_s3(self):
        dataset_params = BogusDatasetParams.get()
        dataset = Dataset.create(self.session, **dataset_params)
        artifact = self.generate_artifact(self.session, dataset.id, upload=True)
        art_bucket_path = artifact.get_bucket_path()
        self.assertS3FileExists(self.bucket, art_bucket_path)
        dataset.asset_deletion()
        self.assertEqual(len(dataset.artifacts), 0)
        self.assertS3FileDoesNotExist(self.bucket, art_bucket_path)

    def test__deployment_directories_delete(self):
        """S3 resources are not deleted"""
        dataset = Dataset.create(self.session, **BogusDatasetParams.get())
        dep_dir = self.generate_deployment_directory(self.session, dataset_id=dataset.id, upload=True)
        cxg_bucket_path = f"{get_cxg_bucket_path(dep_dir)}.cxg/"
        self.assertS3FileExists(self.cellxgene_bucket, cxg_bucket_path)

        dataset.deployment_directories_deletion()
        self.assertS3FileDoesNotExist(self.cellxgene_bucket, cxg_bucket_path)

    def test__dataset_delete(self):
        dataset = Dataset.create(self.session, **BogusDatasetParams.get())
        dataset_id = dataset.id
        artifact = self.generate_artifact(self.session, dataset_id, upload=True)
        art_bucket_path = artifact.get_bucket_path()
        artifact_id = dataset.artifacts[0].id
        dep_dir = self.generate_deployment_directory(self.session, dataset_id=dataset_id, upload=True)
        cxg_bucket_path = f"{get_cxg_bucket_path(dep_dir)}.cxg/"
        deployment_directory_id = dataset.deployment_directories[0].id

        self.session.expire_all()
        dataset = Dataset.get(self.session, dataset_id)
        self.assertIsNotNone(self.session.query(DbDeploymentDirectory).get(deployment_directory_id))
        self.assertIsNotNone(self.session.query(DbDatasetArtifact).get(artifact_id))
        self.assertIsNotNone(dataset)
        self.assertS3FileExists(self.cellxgene_bucket, cxg_bucket_path)
        self.assertS3FileExists(self.bucket, art_bucket_path)
        dataset.delete()
        self.session.expire_all()
        dataset = Dataset.get(self.session, dataset_id)
        self.assertIsNone(dataset)
        self.assertIsNone(self.session.query(DbDeploymentDirectory).get(deployment_directory_id))
        self.assertIsNone(self.session.query(DbDatasetArtifact).get(artifact_id))
        self.assertS3FileExists(self.cellxgene_bucket, cxg_bucket_path)
        self.assertS3FileExists(self.bucket, art_bucket_path)

    def test__generate_tidy_csv_for_all_linked_genesets__correctly_creates_csv(self):
        collection = self.generate_collection(self.session)
        deployment_directory_params = dict(url="test_url")
        dataset = self.generate_dataset(
            self.session, collection=collection, deployment_directories=[deployment_directory_params]
        )
        genes0 = [
            {"gene_symbol": "HBB", "gene_description": "gene 1", "additional_params": {}},
            {"gene_symbol": "HBA1", "gene_description": "gene 2", "additional_params": {}},
            {
                "gene_symbol": "HBA2",
                "gene_description": "gene 3",
                "additional_params": {"provenance1": "some words", "provenance1_description": "another set of words"},
            },
            {"gene_symbol": "HBD", "gene_description": "gene 4", "additional_params": {}},
        ]
        genes1 = [
            {
                "gene_symbol": "IGHG4",
                "gene_description": "gene 1",
                "additional_params": {"provenance1": "some words", "provenance1_description": "another set of words"},
            },
            {
                "gene_symbol": "CANX",
                "gene_description": "gene 2",
                "additional_params": {
                    "provenance1": "some words",
                    "provenance1_description": "another set of words",
                    "provenance2": "some words(2)",
                    "provenance2_description": "another set of words(2)",
                },
            },
            {"gene_symbol": "HBA2", "gene_description": "gene 3", "additional_params": {}},
            {"gene_symbol": "HBD", "gene_description": "gene 4", "additional_params": {}},
        ]
        self.generate_geneset(
            self.session,
            collection=collection,
            dataset_ids=[dataset.id],
            name="first geneset",
            description="describe the geneset",
            genes=genes0,
        )
        self.generate_geneset(
            self.session,
            collection=collection,
            dataset_ids=[dataset.id],
            name="second geneset",
            description="describe another geneset",
            genes=genes1,
        )
        with tempfile.TemporaryDirectory() as temp_dir_name:
            csv_file = dataset.generate_tidy_csv_for_all_linked_genesets(temp_dir_name)
            expected_csv = os.path.abspath(
                os.path.join(os.path.dirname(__file__), "../../../fixtures/sample_geneset_csv.csv")
            )  # noqa
            self.assertTrue(filecmp.cmp(csv_file, expected_csv, shallow=False))

    def test__generate_tidy_csv_for_all_linked_genesets__stores_file_in_correct_location(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        deployment_directory_params = dict(url=f"http://test_domain.example/e/{dataset.id}.cxg/")
        dataset.update(deployment_directories=[deployment_directory_params])
        genes0 = [
            {"gene_symbol": "HBB", "gene_description": "gene 1", "additional_params": {}},
            {"gene_symbol": "HBA1", "gene_description": "gene 2", "additional_params": {}},
            {
                "gene_symbol": "HBA2",
                "gene_description": "gene 3",
                "additional_params": {"provenance1": "some words", "provenance1_description": "another set of words"},
            },
            {"gene_symbol": "HBD", "gene_description": "gene 4", "additional_params": {}},
        ]
        self.generate_geneset(
            self.session,
            collection=collection,
            dataset_ids=[dataset.id],
            name="first geneset",
            description="describe the geneset",
            genes=genes0,
        )
        with tempfile.TemporaryDirectory() as temp_dir_name:
            csv_file = dataset.generate_tidy_csv_for_all_linked_genesets(temp_dir_name)
            s3_file = dataset.copy_csv_to_s3(csv_file)
        expected_suffix = f"{dataset.id}-genesets.csv"
        self.assertTrue(s3_file.endswith(expected_suffix), msg=f"{s3_file} does not end with {expected_suffix}")
        stored_files = [x.key for x in self.cellxgene_bucket.objects.all()]
        self.assertIn(s3_file, stored_files)

        # Delete all geneset files.
        dataset.deployment_directories_deletion()
        self.assertFalse([x.key for x in self.cellxgene_bucket.objects.all()])

    def assertRowsDeleted(self, tests: typing.List[typing.Tuple[str, Base]]):
        """
        Verify if rows have been deleted from the database.
        :param tests: a list of tuples with (primary_key, table)
        """
        self.session.expire_all()
        for p_key, table in tests:
            if len(p_key) == 2:
                # handle the special case for collections with a composite primary key
                actual = [
                    i for i in self.session.query(table).filter(table.id == p_key[0], table.visibility == p_key[1])
                ]
            else:
                actual = [i for i in self.session.query(table).filter(table.id == p_key)]
            self.assertFalse(actual, f"Row not deleted {table.__name__}:{p_key}")

    def assertRowsExist(self, tests):
        """
        Verify if rows exist in the database.
        :param tests: a list of tuples with (primary_key, table)
        """
        self.session.expire_all()
        for p_key, table in tests:
            if len(p_key) == 2:
                # handle the special case for collections with a composite primary key
                actual = [
                    i for i in self.session.query(table).filter(table.id == p_key[0], table.visibility == p_key[1])
                ]
            else:
                actual = [i for i in self.session.query(table).filter(table.id == p_key)]
            self.assertTrue(actual, f"Row does not exist {table.__name__}:{p_key}")
