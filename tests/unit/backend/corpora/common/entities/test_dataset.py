import unittest


from backend.corpora.common.entities.dataset import Dataset


class TestDataset(unittest.TestCase):
    def setUp(self):
        self.uuid = "test_dataset_id"

    def test__get__ok(self):
        dataset = Dataset.get(self.uuid)
        self.assertEqual(dataset.id, self.uuid)
        self.assertEqual(dataset.assay, "test_assay")
        self.assertEqual(dataset.artifacts[0].id, "test_dataset_artifact_id")
        self.assertEqual(dataset.deployment_directories[0].id, "test_deployment_directory_id")
        self.assertEqual(dataset.contributors[0].id, "test_contributor_id")

    def test__get__does_not_exist(self):
        non_existent_key = "non_existent_id"
        self.assertEqual(Dataset.get(non_existent_key), None)
