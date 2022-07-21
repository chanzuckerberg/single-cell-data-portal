from sqlalchemy.orm import Session

from backend.common.entities import Collection, Dataset
from backend.common.entities.geneset import Geneset
from backend.common.entities.collection_link import CollectionLink
from backend.common.utils.db_session import db_session_manager
from tests.unit.backend.utils import (
    BogusCollectionParams,
    BogusDatasetParams,
    BogusGenesetParams,
    BogusDbCollectionLinkParams,
)


class GenerateDataMixin:
    """
    Use to populate the database with test data that should be cleanup after the test
    """

    @staticmethod
    def delete_collection(_id):
        with db_session_manager() as session:
            col = Collection.get(session, _id)
            if col:
                col.delete()

    def generate_collection(self, session: Session, **params) -> Collection:
        _collection = Collection.create(session, **BogusCollectionParams.get(**params))
        self.addCleanup(self.delete_collection, _collection.id)
        return _collection

    @staticmethod
    def delete_dataset(_id):
        with db_session_manager() as session:
            dat = Dataset.get(session, _id)
            if dat:
                dat.delete()

    def generate_dataset(self, session: Session, **params) -> Dataset:
        _dataset = Dataset.create(session, **BogusDatasetParams.get(**params))
        self.addCleanup(self.delete_dataset, _dataset.id)
        return _dataset

    @staticmethod
    def delete_geneset(_id):
        with db_session_manager() as session:
            geneset = Geneset.get(session, _id)
            if geneset:
                geneset.delete()

    def generate_geneset(self, session: Session, **params) -> Geneset:
        _geneset = Geneset.create(session, **BogusGenesetParams.get(**params))
        self.addCleanup(self.delete_geneset, _geneset.id)
        return _geneset

    @staticmethod
    def delete_link(_id):
        with db_session_manager() as session:
            link = CollectionLink.get(session, _id)
            if link:
                link.delete()

    def generate_link(self, session: Session, **params) -> CollectionLink:
        _link = CollectionLink.create(session, **BogusDbCollectionLinkParams.get(**params))
        self.addCleanup(self.delete_link, _link.id)
        return _link
