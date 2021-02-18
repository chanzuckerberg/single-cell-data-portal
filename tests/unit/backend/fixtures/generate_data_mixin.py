from backend.corpora.common.entities import Collection, Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from tests.unit.backend.utils import BogusCollectionParams, BogusDatasetParams


class GenerateDataMixin:
    """
    Use to populate the database with test data that should be cleanup after the test
    """

    def generate_collection(self, session, **params) -> Collection:
        def delete(uuid, visibility):
            with db_session_manager() as session:
                col = Collection.get(session, (uuid, visibility))
                if col:
                    col.delete()

        _collection = Collection.create(session, **BogusCollectionParams.get(**params))
        # Cleanup collection after test
        self.addCleanup(delete, _collection.id, _collection.visibility)
        return _collection

    def generate_dataset(self, session, **params) -> Dataset:
        def delete(uuid):
            with db_session_manager() as session:
                dat = Dataset.get(session, uuid)
                if dat:
                    dat.delete()

        _dataset = Dataset.create(session, **BogusDatasetParams.get(**params))
        # Cleanup collection after test
        self.addCleanup(delete, _dataset.id)
        return _dataset

    def generate_geneset(self, session, **params) -> GeneSet:
        pass
