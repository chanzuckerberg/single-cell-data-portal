from sqlalchemy.exc import IntegrityError

from backend.corpora.common.corpora_orm import DbCollection
from backend.corpora.common.entities import Dataset
from backend.corpora.common.entities.geneset import Geneset, GenesetDatasetLink
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestGeneSets(DataPortalTestCase):
    def test_create_a_geneset(self):
        collection = self.generate_collection(self.session)
        geneset = Geneset.create(
            self.session,
            name="geneset",
            description="words to describe it",
            gene_symbols=["AAA", "BBB", "CCC"],
            collection=collection,
        )
        self.assertEqual(geneset.name, "geneset")

    def test_cant_create_genesets_with_same_name_for_a_collection(self):
        collection = self.generate_collection(self.session)
        self.generate_geneset(session=self.session, name="name", collection=collection)
        with self.assertRaises(IntegrityError):
            self.generate_geneset(session=self.session, name="name", collection=collection)

    def test_retrieve_a_geneset(self):
        geneset = self.generate_geneset(self.session)

        stored_geneset = Geneset.get(self.session, geneset.id)
        self.assertIsInstance(stored_geneset.description, str)
        self.assertIsInstance(stored_geneset.name, str)
        self.assertIsInstance(stored_geneset.gene_symbols, list)
        self.assertGreater(len(stored_geneset.gene_symbols), 1)
        self.assertIsNotNone(stored_geneset.collection)
        self.assertIsInstance(stored_geneset.collection, DbCollection)

    def test_try_to_retrieve_a_geneset_that_does_not_exist(self):
        geneset = Geneset.get(self.session, "not_a_uuid")
        self.assertIsNone(geneset)

    def test_retrieve_all_genesets_for_a_collection(self):
        collection = self.generate_collection(self.session)
        geneset_0 = self.generate_geneset(self.session, collection=collection)
        geneset_1 = self.generate_geneset(self.session, collection=collection)
        geneset_2 = self.generate_geneset(self.session, collection=collection)

        genesets = Geneset.retrieve_all_genesets_for_a_collection(session=self.session, collection_id=collection.id)
        linked_geneset_ids = [x["id"] for x in genesets]
        self.assertIn(geneset_0.id, linked_geneset_ids)
        self.assertIn(geneset_1.id, linked_geneset_ids)
        self.assertIn(geneset_2.id, linked_geneset_ids)


class TestGenesetDatasetLinks(DataPortalTestCase):
    def test_link_a_geneset_and_dataset(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        geneset = self.generate_geneset(session=self.session, collection=collection)
        self.assertEqual(len(geneset.datasets), 0)

        GenesetDatasetLink.create(self.session, dataset_id=dataset.id, geneset_id=geneset.id)
        geneset = Geneset.get(self.session, geneset.id)
        self.assertEqual(len(geneset.datasets), 1)

    def test_get_list_of_genesets_for_dataset(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        self.assertEqual(len(dataset.genesets), 3)

    def test_get_list_of_datasets_for_a_geneset(self):
        collection = self.generate_collection(self.session)
        dataset0 = self.generate_dataset(self.session, collection=collection)
        dataset1 = self.generate_dataset(self.session, collection=collection)
        dataset2 = self.generate_dataset(self.session, collection=collection)
        geneset = self.generate_geneset(
            self.session, collection=collection, dataset_ids=[dataset0.id, dataset1.id, dataset2.id]
        )
        self.assertEqual(len(geneset.datasets), 3)

    def test_deleting_a_geneset_does_not_cascade_to_the_dataset(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        geneset = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        self.assertIsNotNone(dataset.genesets)
        geneset.delete()
        dataset = Dataset.get(self.session, dataset.id)
        self.assertIsNotNone(dataset)
        self.assertEqual(len(dataset.genesets), 0)

    def test_deleting_a_dataset_does_not_cascade_to_the_geneset(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        geneset = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        self.assertIsNotNone(geneset.datasets)
        dataset.delete()
        geneset = Geneset.get(self.session, geneset.id)
        self.assertIsNotNone(geneset)
        self.assertEqual(len(geneset.datasets), 0)
