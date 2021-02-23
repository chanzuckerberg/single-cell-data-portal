from backend.corpora.common.corpora_orm import CollectionVisibility, DbCollection
from backend.corpora.common.entities import Collection
from backend.corpora.common.entities.geneset import Geneset
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestGeneSets(DataPortalTestCase):
    def setUp(self):
        super().setUp()

    def test_create_a_geneset(self):
        test_uuid = "test_genest_id"
        geneset = self.generate_geneset(self.session, id=test_uuid)
        self.assertEqual(test_uuid, geneset.id)

    def test_cant_create_genesets_with_same_name_for_a_collection(self):
        pass

    def test_retrieve_a_geneset(self):
        geneset = self.generate_geneset(self.session)

        stored_geneset = Geneset.get(self.session, geneset.id)
        self.assertIsInstance(stored_geneset.description, str)
        self.assertIsInstance(stored_geneset.name, str)
        self.assertIsInstance(stored_geneset.gene_symbols, list)
        self.assertGreater(len(stored_geneset.gene_symbols), 1)
        self.assertIsNotNone(stored_geneset.collection)
        self.assertIsInstance(stored_geneset.collection, DbCollection)

    def test_retrieve_a_geneset_that_does_not_exist(self):
        geneset = Geneset.get(self.session, "not_a_uuid")
        self.assertIsNone(geneset)

