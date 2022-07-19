from backend.corpora.common.entities.tiledb_data import TileDBData, Utils
from backend.corpora.common.entities.collection_link import CollectionLink
from tests.unit.backend.utils import (
    BogusCollectionParams,
    BogusDatasetParams,
    BogusDbCollectionLinkParams,
)

location = "./test_tiledb/metadata"  # TODO: config this somewhere

class GenerateDataMixin:
    """
    Use to populate the database with test data that should be cleanup after the test
    """

    attrs = Utils.attrs

    @staticmethod
    def delete_collection(_id):
        db = TileDBData(location)
        db.delete_collection(_id)

    def generate_collection(self, **params):
        db = TileDBData(location)
        data = dict(Utils.empty_collection)
        data.update(**BogusCollectionParams.get(**params))
        id = db.create_collection(**data)
        metadata = db.get_collection(id)
        self.addCleanup(self.delete_collection, id)
        return id, metadata

    @staticmethod
    def delete_dataset(coll_id, dataset_id):
        db = TileDBData(location)
        db.delete_dataset(coll_id, dataset_id)

    def generate_dataset(self, coll_id: str, **params):
        db = TileDBData(location)
        data = dict(Utils.empty_dataset)
        data.update(**BogusDatasetParams.get(**params))
        dataset_id = db.add_dataset(coll_id, data)
        dataset = db.get_dataset(dataset_id)
        self.addCleanup(self.delete_dataset, dataset_id)
        return dataset_id, dataset

    # @staticmethod
    # def delete_geneset(_id):
    #     with db_session_manager() as session:
    #         geneset = Geneset.get(session, _id)
    #         if geneset:
    #             geneset.delete()

    # def generate_geneset(self, session: Session, **params) -> Geneset:
    #     _geneset = Geneset.create(session, **BogusGenesetParams.get(**params))
    #     self.addCleanup(self.delete_geneset, _geneset.id)
    #     return _geneset

    def generate_link(self, session: Session, **params) -> CollectionLink:
        _link = BogusDbCollectionLinkParams.get(**params)
        return _link
