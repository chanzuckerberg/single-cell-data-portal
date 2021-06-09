import filecmp
import os
import tempfile

from tests.unit.backend.corpora.common.entities.datasets import TestDataset


class TestDatasetGenesetCSV(TestDataset):
    def test__generate_tidy_csv_for_all_linked_genesets__correctly_creates_csv(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        genes0 = [
            {"gene_symbol": "HBB", "gene_description": "gene 1", "additional_params": {}},
            {"gene_symbol": "HBA1", "gene_description": "gene 2", "additional_params": {}},
            {
                "gene_symbol": "HBA2",
                "gene_description": "gene 3",
                "additional_params": {"provenance1": "some words", "provenance1_description": "another set of words"},
            },
            {"gene_symbol": "HBD", "gene_description": "gene 4", "additional_params": {}},
        ]
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
        self.generate_geneset(
            self.session,
            collection=collection,
            dataset_ids=[dataset.id],
            name="first geneset",
            description="describe the geneset",
            genes=genes0,
        )
        self.generate_geneset(
            self.session,
            collection=collection,
            dataset_ids=[dataset.id],
            name="second geneset",
            description="describe another geneset",
            genes=genes1,
        )
        with tempfile.TemporaryDirectory() as temp_dir_name:
            csv_file = dataset.generate_tidy_csv_for_all_linked_genesets(temp_dir_name)
            expected_csv = os.path.abspath(
                os.path.join(os.path.dirname(__file__), "../../../../fixtures/sample_geneset_csv.csv")
            )  # noqa
            self.assertTrue(filecmp.cmp(csv_file, expected_csv, shallow=False))

    def test__generate_tidy_csv_for_all_linked_genesets__stores_file_in_correct_location(self):
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset(self.session, collection=collection)
        dataset.update(explorer_url=f"http://test_domain.example/e/{dataset.id}.cxg/")
        genes0 = [
            {"gene_symbol": "HBB", "gene_description": "gene 1", "additional_params": {}},
            {"gene_symbol": "HBA1", "gene_description": "gene 2", "additional_params": {}},
            {
                "gene_symbol": "HBA2",
                "gene_description": "gene 3",
                "additional_params": {"provenance1": "some words", "provenance1_description": "another set of words"},
            },
            {"gene_symbol": "HBD", "gene_description": "gene 4", "additional_params": {}},
        ]
        self.generate_geneset(
            self.session,
            collection=collection,
            dataset_ids=[dataset.id],
            name="first geneset",
            description="describe the geneset",
            genes=genes0,
        )
        with tempfile.TemporaryDirectory() as temp_dir_name:
            csv_file = dataset.generate_tidy_csv_for_all_linked_genesets(temp_dir_name)
            s3_file = dataset.copy_csv_to_s3(csv_file)
        expected_suffix = f"{dataset.id}-genesets.csv"

        self.assertTrue(s3_file.endswith(expected_suffix), msg=f"{s3_file} does not end with {expected_suffix}")
        stored_files = [x.key for x in self.cellxgene_bucket.objects.all()]
        self.assertIn(s3_file, stored_files)

        # Delete all geneset files.

        dataset.delete_explorer_cxg_object_from_s3()
        self.assertFalse([x.key for x in self.cellxgene_bucket.objects.all()])
