from time import sleep

from backend.corpora.common.entities import Dataset
from tests.unit.backend.corpora.common.entities.datasets import TestDataset


class TestGetDatasetByExplorerUrl(TestDataset):
    def test__get__ok(self):
        created_dataset = self.generate_dataset(self.session)
        dataset = Dataset.get_by_explorer_url(self.session, created_dataset.explorer_url)
        self.assertEqual(created_dataset.id, dataset.id)

    def test__get__ok_with_final_slash(self):
        created_dataset = self.generate_dataset(self.session, explorer_url="test_url2")
        dataset = Dataset.get_by_explorer_url(self.session, "test_url2/")
        self.assertEqual(created_dataset.id, dataset.id)

    def test__get__ok_without_final_slash(self):
        created_dataset = self.generate_dataset(self.session, explorer_url="test_url3/")
        dataset = Dataset.get_by_explorer_url(self.session, "test_url3")
        self.assertEqual(created_dataset.id, dataset.id)

    def test__get__ok_with_consistent_slashes(self):
        created_dataset = self.generate_dataset(self.session, explorer_url="test_url4/")
        dataset = Dataset.get_by_explorer_url(self.session, "test_url4/")
        self.assertEqual(created_dataset.id, dataset.id)

    def test__get__does_not_exist(self):
        non_existent_url = "made_up_url"
        self.assertIsNone(Dataset.get_by_explorer_url(self.session, non_existent_url))

    def test__get__url_associated_with_most_recently_created_datasets(self):
        url = "some_url"
        self.generate_dataset(self.session, name="older", explorer_url=url)
        sleep(0.01)
        dataset_1 = self.generate_dataset(self.session, name="younger", explorer_url=url)
        dataset = Dataset.get_by_explorer_url(self.session, url)
        self.assertEqual(dataset.name, "younger")
        self.assertEqual(dataset.id, dataset_1.id)
