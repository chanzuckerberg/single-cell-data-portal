from backend.corpora.common.entities import Collection, Dataset
from tests.unit.backend.utils import BogusCollectionParams, BogusDatasetParams


class GenerateDataMixin:
    """
    Use to populate the database with test data that should be cleanup after the test
    """

    def generate_collection(self, **params) -> Collection:
        _collection = Collection.create(**BogusCollectionParams.get(**params))
        # Cleanup collection after test
        self.addCleanup(_collection.delete)
        return _collection

    def generate_dataset(self, **params) -> Dataset:
        _dataset = Dataset.create(**BogusDatasetParams.get(**params))
        # Cleanup collection after test
        self.addCleanup(_dataset.delete)
        return _dataset
