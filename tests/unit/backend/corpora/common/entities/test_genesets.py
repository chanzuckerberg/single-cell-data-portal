from sqlalchemy.exc import IntegrityError

from backend.corpora.common.corpora_orm import DbCollection, generate_uuid
from backend.corpora.common.entities import Dataset
from backend.corpora.common.entities.geneset import Geneset, GenesetDatasetLink
from backend.corpora.common.utils.exceptions import CorporaException
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestGeneSets(DataPortalTestCase):
    def test_create_a_geneset(self):
        collection = self.generate_collection(self.session)
        geneset = Geneset.create(
            self.session,
            name="geneset",
            description="words to describe it",
            genes=["AAA", "BBB", "CCC"],
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
        self.assertIsInstance(stored_geneset.genes, list)
        self.assertGreater(len(stored_geneset.genes), 1)
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

    def test_gene_order_is_maintained(self):
        collection = self.generate_collection(self.session)
        genes = [
            {"name": "a", "description": "words"},
            {"name": "b", "description": "words"},
            {"name": "c", "description": "words"},
            {"name": "d", "description": "words"},
        ]
        geneset = self.generate_geneset(session=self.session, name="name", genes=genes, collection=collection)
        self.session.flush()
        geneset = Geneset.get(self.session, geneset.id)
        self.assertEqual(geneset.genes[0], genes[0])
        self.assertEqual(geneset.genes, genes)

    def test_geneset_to_gene_dict_conversion__ok(self):
        expected_gene_dicts = [
            {
                "GENE_SET_NAME": "first geneset",
                "GENE_SET_DESCRIPTION": "describe the geneset",
                "GENE_SYMBOL": "IGHG4",
                "GENE_DESCRIPTION": "gene 1",
                "PROVENANCE1": "some words",
                "PROVENANCE1_DESCRIPTION": "another set of words",
            },
            {
                "GENE_SET_NAME": "first geneset",
                "GENE_SET_DESCRIPTION": "describe the geneset",
                "GENE_SYMBOL": "CANX",
                "GENE_DESCRIPTION": "gene 2",
                "PROVENANCE1": "some words",
                "PROVENANCE2": "some words(2)",
                "PROVENANCE1_DESCRIPTION": "another set of words",
                "PROVENANCE2_DESCRIPTION": "another set of words(2)",
            },
            {
                "GENE_SET_NAME": "first geneset",
                "GENE_SET_DESCRIPTION": "describe the geneset",
                "GENE_SYMBOL": "HBA2",
                "GENE_DESCRIPTION": "gene 3",
            },
            {
                "GENE_SET_NAME": "first geneset",
                "GENE_SET_DESCRIPTION": "describe the geneset",
                "GENE_SYMBOL": "HBD",
                "GENE_DESCRIPTION": "gene 4",
            },
        ]
        collection = self.generate_collection(self.session)
        genes1 = [
            {
                "gene_symbol": "IGHG4",
                "gene_description": "gene 1",
                "additional_params": {"provenance1": "some words", "provenance1_description": "another set of words"},
            },
            {
                "gene_symbol": "CANX",
                "gene_description": "gene 2",
                "additional_params": {
                    "provenance1": "some words",
                    "provenance1_description": "another set of words",
                    "provenance2": "some words(2)",
                    "provenance2_description": "another set of words(2)",
                },
            },
            {"gene_symbol": "HBA2", "gene_description": "gene 3", "additional_params": {}},
            {"gene_symbol": "HBD", "gene_description": "gene 4", "additional_params": {}},
        ]
        geneset = self.generate_geneset(
            self.session, collection=collection, name="first geneset", description="describe the geneset", genes=genes1
        )

        gene_dict, max_parms = geneset.convert_geneset_to_gene_dicts()
        with self.subTest("geneset to gene dict conversion returns expected list"):
            self.assertEqual(expected_gene_dicts, gene_dict)

        with self.subTest("geneset to gene dict conversion correctly counts max params"):
            self.assertEqual(max_parms, 2)


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
        geneset_dataset_link = GenesetDatasetLink.get(self.session, geneset_id=geneset.id, dataset_id=dataset.id)
        self.assertIsNotNone(geneset_dataset_link)
        geneset.delete()
        dataset = Dataset.get(self.session, dataset.id)
        self.assertIsNotNone(dataset)
        self.assertEqual(len(dataset.genesets), 0)
        geneset_dataset_link = GenesetDatasetLink.get(self.session, geneset_id=geneset.id, dataset_id=dataset.id)
        self.assertIsNone(geneset_dataset_link)

    def test_deleting_a_dataset_does_not_cascade_to_the_geneset(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        geneset = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        self.assertIsNotNone(geneset.datasets)
        dataset.delete()
        geneset = Geneset.get(self.session, geneset.id)
        self.assertIsNotNone(geneset)
        self.assertEqual(len(geneset.datasets), 0)

    def test_get_gene_set_dataset_link__ok(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        geneset = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        link = GenesetDatasetLink.get(self.session, dataset_id=dataset.id, geneset_id=geneset.id)
        self.assertIsInstance(link, GenesetDatasetLink)

    def test_get_gene_set_dataset_link__does_not_exist(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        geneset = self.generate_geneset(self.session, collection=collection)
        link = GenesetDatasetLink.get(self.session, dataset_id=dataset.id, geneset_id=geneset.id)
        self.assertIsNone(link)

    def test__update_links_for_a_dataset__ok(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        geneset0 = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        geneset1 = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        geneset2 = self.generate_geneset(self.session, collection=collection)
        GenesetDatasetLink.update_links_for_a_dataset(
            self.session, dataset_id=dataset.id, add=[geneset2.id], remove=[geneset0.id]
        )

        linked_genesest_ids = [x.id for x in dataset.genesets]
        self.assertIn(geneset1.id, linked_genesest_ids)
        self.assertIn(geneset2.id, linked_genesest_ids)
        self.assertNotIn(geneset0.id, linked_genesest_ids)

    def test__update_links_for_a_dataset__if_one_fails_no_links_are_changed(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        geneset0 = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        geneset1 = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        geneset2 = self.generate_geneset(self.session, collection=collection)
        with self.subTest("cant link geneset that doesnot exist"):
            with self.assertRaises(CorporaException):
                GenesetDatasetLink.update_links_for_a_dataset(
                    self.session, dataset_id=dataset.id, add=[geneset2.id, generate_uuid()], remove=[geneset0.id]
                )

            linked_genesest_ids = [x.id for x in dataset.genesets]
            self.assertIn(geneset0.id, linked_genesest_ids)
            self.assertIn(geneset1.id, linked_genesest_ids)
            self.assertNotIn(geneset2.id, linked_genesest_ids)
        with self.subTest("cant delete link to geneset that doesnot exist"):
            with self.assertRaises(CorporaException):
                GenesetDatasetLink.update_links_for_a_dataset(
                    self.session, dataset_id=dataset.id, add=[geneset2.id], remove=[geneset0.id, generate_uuid()]
                )

            linked_genesest_ids = [x.id for x in dataset.genesets]
            self.assertIn(geneset0.id, linked_genesest_ids)
            self.assertIn(geneset1.id, linked_genesest_ids)
            self.assertNotIn(geneset2.id, linked_genesest_ids)
