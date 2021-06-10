from backend.corpora.common.entities import Dataset
from tests.unit.backend.corpora.common.entities.datasets import TestDataset


class TestGetDataset(TestDataset):
    def test__get__ok(self):

        self.create_dataset_with_artifacts()
        dataset = Dataset.get(self.session, self.uuid)
