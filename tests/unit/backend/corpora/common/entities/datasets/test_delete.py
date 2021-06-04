from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
    DbDataset,
    DbDatasetArtifact,
    DbCollection,
    UploadStatus,
)
from backend.corpora.common.entities import Dataset
from backend.corpora.common.entities.dataset import get_cxg_bucket_path
from backend.corpora.common.entities.geneset import GenesetDatasetLink, Geneset
from tests.unit.backend.corpora.common.entities.datasets import TestDataset
from tests.unit.backend.utils import BogusDatasetParams


class TestDeleteDataset(TestDataset):
    def test__cascade_delete_dataset__ok(self):
        # Create the dataset
        test_dataset = Dataset.create(
            self.session,
            **BogusDatasetParams.get(
                collection_id="test_collection_id",
                collection_visibility=CollectionVisibility.PUBLIC.name,
                artifacts=[{}],
            ),
        )
        test_dataset_ids = [(test_dataset.id, DbDataset)]
        test_artifact_ids = [(art.id, DbDatasetArtifact) for art in test_dataset.artifacts]
        test_collection_ids = [(("test_collection_id", CollectionVisibility.PUBLIC.name), DbCollection)]

        with self.subTest("verify everything exists"):
            expected_exists = test_collection_ids + test_dataset_ids + test_artifact_ids
            self.assertRowsExist(expected_exists)

        # Delete the dataset
        test_dataset.delete()

        with self.subTest("Verify Deletion"):
            expected_deleted = test_dataset_ids + test_artifact_ids
            expected_exists = test_collection_ids
            self.assertRowsDeleted(expected_deleted)
            self.assertRowsExist(expected_exists)

    def test__tombstone_dataset_and_delete_child_objects(self):
        dataset = self.create_dataset_with_artifacts(artifact_count=3)
        geneset = self.generate_geneset(self.session)
        GenesetDatasetLink.create(self.session, geneset.id, dataset.id)
        self.assertEqual(dataset.processing_status.upload_status, UploadStatus.UPLOADING)
        self.assertEqual(len(dataset.artifacts), 3)
        self.assertEqual(len(dataset.genesets), 1)
        self.assertFalse(dataset.tombstone)

        dataset.tombstone_dataset_and_delete_child_objects()
        self.assertEqual(len(dataset.artifacts), 0)
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

    def test__delete_explorer_cxg_object_from_s3(self):
        dataset = Dataset.create(self.session, **BogusDatasetParams.get())
        self.create_explorer_s3_object(self.session, dataset_id=dataset.id, upload=True)
        cxg_bucket_path = f"{get_cxg_bucket_path(dataset.explorer_url)}.cxg/"
        self.assertS3FileExists(self.cellxgene_bucket, cxg_bucket_path)

        dataset.delete_explorer_cxg_object_from_s3()
        self.assertS3FileDoesNotExist(self.cellxgene_bucket, cxg_bucket_path)

    def test__dataset_delete(self):
        dataset = Dataset.create(self.session, **BogusDatasetParams.get())
        dataset_id = dataset.id
        artifact = self.generate_artifact(self.session, dataset_id, upload=True)
        art_bucket_path = artifact.get_bucket_path()
        artifact_id = dataset.artifacts[0].id
        self.create_explorer_s3_object(self.session, dataset_id=dataset.id, upload=True)
        cxg_bucket_path = f"{get_cxg_bucket_path(dataset.explorer_url)}.cxg/"

        self.session.expire_all()
        dataset = Dataset.get(self.session, dataset_id)
        self.assertIsNotNone(dataset.explorer_url)
        self.assertIsNotNone(self.session.query(DbDatasetArtifact).get(artifact_id))
        self.assertIsNotNone(dataset)
        self.assertS3FileExists(self.cellxgene_bucket, cxg_bucket_path)
        self.assertS3FileExists(self.bucket, art_bucket_path)
        dataset.delete()
        self.session.expire_all()
        dataset = Dataset.get(self.session, dataset_id)
        self.assertIsNone(dataset)
        self.assertIsNone(self.session.query(DbDatasetArtifact).get(artifact_id))
        self.assertS3FileExists(self.cellxgene_bucket, cxg_bucket_path)
        self.assertS3FileExists(self.bucket, art_bucket_path)
