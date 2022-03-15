import os
import shutil
import unittest
from unittest.mock import patch

import anndata
import tiledb

from backend.wmg.data.cube_pipeline import load
from backend.wmg.data.load_corpus import load_h5ad
from backend.wmg.data.schemas.corpus_schema import create_tdb
from tests.unit.backend.wmg.fixtures.test_anndata_object import create_anndata_test_object
relative_fixtures_path = f"../fixtures/"


class TestCorpusLoad(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        super().setUp(cls)
        cls.tmp_file = f"{relative_fixtures_path}tmp"
        os.mkdir(cls.tmp_file)

        basic_test_anndata_object = create_anndata_test_object(num_genes=3, num_cells=5)
        larger_test_anndata_object = create_anndata_test_object(num_genes=1000, num_cells=5000)

        cls.path_to_basic_anndata_test_object = f"{cls.tmp_file}/basic_test_dataset.h5ad"
        cls.path_to_larger_test_anndata_object = f"{cls.tmp_file}/larger_test_dataset.h5ad"

        basic_test_anndata_object.write(cls.path_to_basic_anndata_test_object, compression="gzip")
        larger_test_anndata_object.write(cls.path_to_larger_test_anndata_object, compression="gzip")

    @classmethod
    def tearDownClass(cls) -> None:
        super().tearDownClass()
        shutil.rmtree(cls.tmp_file)

    def setUp(self) -> None:
        super().setUp()
        self.path_to_dataset_0 = f"{relative_fixtures_path}small_datasets/be573215-4116-40a1-9c0f-1b21ba53482b/local.h5ad"
        self.path_to_dataset_1 = f"{relative_fixtures_path}small_datasets/d50b8959-6ce9-4a9b-b804-99892c93b183/local.h5ad"
        self.path_to_datasets = f"{relative_fixtures_path}small_datasets"
        self.corpus_name = f"{self.tmp_file}/test-group"

        if not tiledb.VFS().is_dir(self.corpus_name):
            create_tdb(self.corpus_name)
        # self.pbmc3k_anndata_object = get_test_anndata_dataset()

    def tearDown(self) -> None:
        super().tearDown()
        shutil.rmtree(self.corpus_name)

    @patch('backend.wmg.data.cube_pipeline.tiledb.consolidate')
    @patch('backend.wmg.data.cube_pipeline.tiledb.vacuum')
    @patch("backend.wmg.data.cube_pipeline.load_h5ad")
    def test__load_loads_all_datasets_in_directory(self, mock_load_h5ad, mock_vacuum, mock_consolidate):
        load(self.path_to_datasets, self.corpus_name)
        self.assertEqual(mock_load_h5ad.call_count, 2)
        self.assertEqual(mock_vacuum.call_count, 4)
        self.assertEqual(mock_consolidate.call_count, 4)

    @patch('backend.wmg.data.load_corpus.update_global_var')
    @patch('backend.wmg.data.load_corpus.validate_dataset_properties')
    def test_invalid_datasets_are_not_added_to_corpus(self, mock_validation, mock_global_var):
        mock_validation.return_value = False
        load(self.path_to_datasets, self.corpus_name)
        self.assertEqual(mock_global_var.call_count, 0)

    def test_global_var_array_updated_when_dataset_contains_new_genes(self):
        load_h5ad(self.path_to_dataset_0, self.corpus_name, False)
        load_h5ad(self.path_to_dataset_1, self.corpus_name, False)

        with tiledb.open(f"{self.corpus_name}/var", "r") as var:
            var_df = var.df[:]
            total_stored_genes = len(set(var_df['gene_ontology_term_id'].to_numpy(dtype=str)))
        total_genes = [12]
        self.assertEqual(len(total_genes), total_stored_genes)

    def test_mapping_between_local_file_and_global_tdb_is_valid_and_consistent_as_datasets_are_added(self):
        """
        DO NOT DELETE THIS TEST
        """
        # load_h5ad(self.path_to_dataset_0, self.corpus_name, False)
        pass




    def test_axes_labels_updated_for_new_genes(self):
        pass

    def test_axes_labels_saved_as_meta_data_for_tdb_object(self):
        pass

    def test_raw_expression_matrix_saved_as_sparse_matrix(self):
        pass

    def test_expression_matrix_saved_as_sparse_matrix(self):
        pass

    def test_expression_matrix_matches_global_index(self):
        pass

    def test_raw_expression_matrix_normalized_by_rankit(self):
        pass
