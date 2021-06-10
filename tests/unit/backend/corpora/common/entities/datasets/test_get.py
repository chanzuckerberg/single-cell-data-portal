from backend.corpora.common.corpora_orm import DbDatasetArtifact, DbDatasetProcessingStatus
from backend.corpora.common.entities import Dataset
from tests.unit.backend.corpora.common.entities.datasets import TestDataset


class TestGetDataset(TestDataset):
    def test__get__ok(self):
        dataset = Dataset.get(self.session, self.uuid)
        self.assertEqual(dataset.id, self.uuid)
        self.assertEqual(len(dataset.assay), 1)
        self.assertDictEqual(dataset.assay[0], {"ontology_term_id": "test_obo", "label": "test_assay"})

        # Verify Artifact relationship
        self.assertIsInstance(dataset.artifacts[0], DbDatasetArtifact)
        self.assertEqual(dataset.artifacts[0].id, "test_dataset_artifact_id")

        # Verify Processing Status relationship
        self.assertIsInstance(dataset.processing_status, DbDatasetProcessingStatus)
        self.assertEqual(dataset.processing_status.id, "test_dataset_processing_status_id")
        self.assertEqual(dataset.processing_status.dataset_id, "test_dataset_id")

    def test__get__does_not_exist(self):
        non_existent_key = "non_existent_id"
        self.assertEqual(Dataset.get(self.session, non_existent_key), None)

    def test__get_asset__ok(self):
        dataset = Dataset.get(self.session, self.uuid)
        expected_asset_id = "test_dataset_artifact_id"
        asset = dataset.get_asset("test_dataset_artifact_id")
        self.assertEqual(expected_asset_id, asset.id)

    def test__get_asset__not_found(self):
        dataset = Dataset.get(self.session, self.uuid)
        asset = dataset.get_asset("fake_asset")
        self.assertIsNone(asset)
