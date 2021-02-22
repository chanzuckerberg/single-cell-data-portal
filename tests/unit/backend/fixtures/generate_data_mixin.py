import random
import string

from backend.corpora.common.entities import Collection, Dataset
from backend.corpora.common.entities.geneset import Geneset
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

    def generate_geneset(self, session, **params) -> Geneset:
        def delete(uuid):
            with db_session_manager() as session:
                geneset = Geneset.get(session, uuid)
                if geneset:
                    geneset.delete()

        bogus_description = "This is a geneset bwhahaha"
        bogus_name = self.generate_random_string(7)
        gene_symbols = []
        for i in range(6):
            gene_symbols.append(self.generate_random_string(4))
        _geneset = Geneset.create(session, name=bogus_name, description=bogus_description, gene_symbols=gene_symbols)

    @staticmethod
    def generate_random_string(self, length=7):
        return ''.join(random.choice(string.ascii_letters) for i in range(length))
