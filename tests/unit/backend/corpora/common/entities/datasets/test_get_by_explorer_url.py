from backend.corpora.common.entities import Dataset
from tests.unit.backend.corpora.common.entities.datasets import TestDataset


class TestGetDatasetByExplorerUrl(TestDataset):
    def test__get__ok(self):
        created_dataset = self.generate_dataset(self.session)
        dataset = Dataset.get_by_explorer_url(created_dataset.explorer_url)
        self.assertEqual(created_dataset.id, dataset.id)

    def test__get__does_not_exist(self):
        non_existent_url = "made_up_url"
        self.assertIsNone(Dataset.get_by_explorer_url(non_existent_url))

    def test__get__url_associated_with_multiple_datasets(self):
        url = "some_url"
        import pdb
        pdb.set_trace()
        self.generate_dataset(self.session, explorer_url=url)
        self.generate_dataset(self.session, explorer_url=url)
        dataset = Dataset.get_by_explorer_url(url)

        assert 1 == 2
