from backend.corpora.common.entities.tiledb_data import db, Utils
from backend.corpora.common.entities.collection_link import CollectionLink
from tests.unit.backend.utils import (
    BogusCollectionParams,
    BogusDatasetParams,
    BogusDbCollectionLinkParams,
)

class GenerateDataMixin:
    """
    Use to populate the database with test data that should be cleanup after the test
    """

    attrs = Utils.attrs

    def delete_collection(self, _id):
        db.delete_collection(_id)

    def generate_collection(self, **params):
        
        data = dict(Utils.empty_collection)
        data.update(**BogusCollectionParams.get(**params))
        id = db.create_collection(data)
        metadata = db.get_collection(id)
        self.addCleanup(self.delete_collection, id)
        return id, metadata

    def get_collection(self, coll_id: str):
        
        return db.get_collection(coll_id)

    def get_doi(self, links) -> str:
        doi = [link for link in links if link['link_type'] == "DOI"]
        if doi:
            return doi[0]['link_url']
        else:
            return None

    def publish_collection(self, coll_id: str):
        
        db.publish_collection(coll_id)

    def revise_collection(self, coll_id: str):
        
        id = db.create_revision(coll_id)
        return id

    def delete_dataset(self, coll_id, dataset_id):
        
        db.delete_dataset(coll_id, dataset_id)

    def generate_dataset(self, coll_id: str, **params):
        
        data = dict(Utils.empty_dataset)
        data.update(**BogusDatasetParams.get(**params))
        dataset_id = db.add_dataset(coll_id, data)
        dataset = db.get_dataset(dataset_id)
        self.addCleanup(self.delete_dataset, coll_id, dataset_id)
        return dataset_id, dataset
    
    def get_dataset(self, dataset_id: str):
        
        return db.get_dataset(dataset_id)

    def update_processing(self, dataset_id: str, status: dict):
        
        db.edit_dataset(dataset_id, "processing_status", status)

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

    def generate_link(self, **params) -> CollectionLink:
        _link = BogusDbCollectionLinkParams.get(**params)
        return _link
