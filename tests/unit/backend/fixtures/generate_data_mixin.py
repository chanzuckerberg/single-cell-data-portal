from backend.corpora.common.entities import Collection, Dataset
from tests.unit.backend.utils import BogusCollectionParams, BogusDatasetParams


class GenerateDataMixin:
    """
    Use to populate the database with test data that should be cleanup after the test
    """

    def generate_collection(self, **params) -> Collection:
        def delete(uuid, visibility):
            col = Collection.get((uuid, visibility))
            if col:
                col.delete()

        _collection = Collection.create(**BogusCollectionParams.get(**params))
        # Cleanup collection after test
        self.addCleanup(delete, _collection.id, _collection.visibility)
        return _collection

    def generate_dataset(self, **params) -> Dataset:
        def delete(uuid):
            dat = Dataset.get(uuid)
            if dat:
                dat.delete()

        _dataset = Dataset.create(**BogusDatasetParams.get(**params))
        # Cleanup collection after test
        self.addCleanup(delete, _dataset.id)
        return _dataset
