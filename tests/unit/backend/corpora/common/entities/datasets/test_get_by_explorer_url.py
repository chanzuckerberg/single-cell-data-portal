from time import sleep

from backend.corpora.common.entities import Dataset
from tests.unit.backend.corpora.common.entities.datasets import TestDataset


class TestGetDatasetByExplorerUrl(TestDataset):
    def test__get__ok(self):
        created_dataset = self.generate_dataset(self.session)
        dataset = Dataset.get_by_explorer_url(self.session, created_dataset.explorer_url)
        self.assertEqual(created_dataset.id, dataset.id)

    def test__get__does_not_exist(self):
        non_existent_url = "made_up_url"
        self.assertIsNone(Dataset.get_by_explorer_url(self.session, non_existent_url))

    def test__get__url_associated_with_most_recently_created_datasets(self):
        url = "some_url"

        dataset_0 = self.generate_dataset(self.session, name="older", explorer_url=url)
        sleep(.01)
        dataset_1 = self.generate_dataset(self.session, name="younger", explorer_url=url)
        dataset = Dataset.get_by_explorer_url(self.session, url)
        self.assertEqual(dataset.name, "younger")
        self.assertEqual(dataset.id, dataset_1.id)
