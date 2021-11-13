from backend.corpora.common.entities import Collection, Dataset, DatasetAsset
from backend.corpora.common.entities.geneset import Geneset
from backend.corpora.common.utils.db_session import db_session_manager
from tests.unit.backend.utils import (
    BogusCollectionParams,
    BogusDatasetParams,
    BogusGenesetParams,
    BogusDatasetAssetParams,
)


class GenerateDataMixin:
    """
    Use to populate the database with test data that should be cleanup after the test
    """

    @staticmethod
    def delete_collection(uuid, visibility):
        with db_session_manager() as session:
            col = Collection.get(session, (uuid, visibility))
            if col:
                col.delete()

    def generate_collection(self, session, **params) -> Collection:
        _collection = Collection.create(session, **BogusCollectionParams.get(**params))
        self.addCleanup(self.delete_collection, _collection.id, _collection.visibility)
        return _collection

    @staticmethod
    def delete_dataset(uuid):
        with db_session_manager() as session:
            dat = Dataset.get(session, uuid)
            if dat:
                dat.delete()

    def generate_dataset(self, session, **params) -> Dataset:
        _dataset = Dataset.create(session, **BogusDatasetParams.get(**params))
        self.addCleanup(self.delete_dataset, _dataset.id)
        return _dataset

    def generate_asset(self, session, dataset_id, **params) -> Dataset:
        _dataset = DatasetAsset.create(session, dataset_id, **BogusDatasetAssetParams.get(**params))
        return _dataset

    @staticmethod
    def delete_geneset(uuid):
        with db_session_manager() as session:
            geneset = Geneset.get(session, uuid)
            if geneset:
                geneset.delete()

    def generate_geneset(self, session, **params) -> Geneset:
        _geneset = Geneset.create(session, **BogusGenesetParams.get(**params))
        self.addCleanup(self.delete_geneset, _geneset.id)
        return _geneset
