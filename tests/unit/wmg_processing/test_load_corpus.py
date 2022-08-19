import os
import pathlib
import shutil
import tempfile
import unittest
from unittest.mock import patch

import anndata
import numpy as np
import tiledb
from scipy import sparse
from scipy.sparse import coo_matrix, csr_matrix
from backend.wmg.data.rankit import rankit
from backend.wmg.data.cube_pipeline import load_data_and_create_cube
from backend.corpus_asset_pipelines.integrated_corpus.job import build_integrated_corpus
from backend.corpus_asset_pipelines.integrated_corpus.load import load_dataset
from backend.corpus_asset_pipelines.integrated_corpus.validate import validate_dataset_properties
from backend.corpus_asset_pipelines.integrated_corpus.transform import (
    filter_out_rankits_with_low_expression_counts,
    apply_pre_concatenation_filters,
)
from backend.wmg.data.constants import RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD
from backend.wmg.data.schemas.corpus_schema import create_tdb_integrated_corpus, OBS_ARRAY_NAME, VAR_ARRAY_NAME
from tests.unit.backend.wmg.fixtures.test_anndata_object import create_anndata_test_object


class TestCorpusLoad(unittest.TestCase):
    @staticmethod
    def fixture_file_path(relative_filename):
        import os

        parent_dir_of_curr_file = "/".join(os.path.dirname(__file__).split("/")[:-1])
        return os.path.abspath(os.path.join(parent_dir_of_curr_file, relative_filename))

    @classmethod
    def setUpClass(cls) -> None:
        super().setUp(cls)
        cls.tmp_dir = tempfile.mkdtemp()

        basic_test_anndata_object = create_anndata_test_object(num_genes=3, num_cells=5)
        larger_test_anndata_object = create_anndata_test_object(num_genes=1000, num_cells=5000)
        os.mkdir(f"{cls.tmp_dir}/datasets")
        os.mkdir(f"{cls.tmp_dir}/datasets/basic_test_dataset")
        os.mkdir(f"{cls.tmp_dir}/datasets/larger_test_dataset")

        cls.small_anndata_filename = pathlib.Path(cls.tmp_dir, "datasets/basic_test_dataset/local.h5ad")
        cls.large_anndata_filename = pathlib.Path(cls.tmp_dir, "datasets/larger_test_dataset/local.h5ad")

        cls.small_anndata_filename.touch()
        cls.large_anndata_filename.touch()

        basic_test_anndata_object.write(cls.small_anndata_filename, compression="gzip")
        larger_test_anndata_object.write(cls.large_anndata_filename, compression="gzip")

    @classmethod
    def tearDownClass(cls) -> None:
        super().tearDownClass()
        shutil.rmtree(cls.tmp_dir)

    def setUp(self) -> None:
        super().setUp()
        self.path_to_datasets = pathlib.Path(self.tmp_dir, "datasets")
        self.corpus_name = "test-group"
        self.corpus_path = f"{self.tmp_dir}/{self.corpus_name}"
        if not tiledb.VFS().is_dir(self.corpus_path):
            create_tdb_integrated_corpus(self.corpus_path)
            # self.pbmc3k_anndata_object = get_test_anndata_dataset()

    def tearDown(self) -> None:
        super().tearDown()
        shutil.rmtree(self.corpus_path)

    @patch("backend.corpus_asset_pipelines.integrated_corpus.job.tiledb.consolidate")
    @patch("backend.corpus_asset_pipelines.integrated_corpus.job.tiledb.vacuum")
    @patch("backend.corpus_asset_pipelines.integrated_corpus.job.process_h5ad_for_corpus")
    def test__load_loads_all_datasets_in_directory(self, mock_process_h5ad, mock_vacuum, mock_consolidate):
        build_integrated_corpus(self.path_to_datasets, self.corpus_path)
        self.assertEqual(mock_process_h5ad.call_count, 2)
        self.assertEqual(mock_vacuum.call_count, 3)
        self.assertEqual(mock_consolidate.call_count, 3)

    @patch("backend.corpus_asset_pipelines.integrated_corpus.load.update_corpus_var")
    @patch("backend.corpus_asset_pipelines.integrated_corpus.job.validate_dataset_properties")
    def test_invalid_datasets_are_not_added_to_corpus(self, mock_validation, mock_global_var):
        mock_validation.return_value = False
        build_integrated_corpus(self.path_to_datasets, self.corpus_path)
        self.assertEqual(mock_global_var.call_count, 0)

    @patch("backend.corpus_asset_pipelines.integrated_corpus.load.transform_dataset_raw_counts_to_rankit")
    def test_global_var_array_updated_when_dataset_contains_new_genes(self, mock_rankit):
        small_anndata_object = anndata.read_h5ad(self.small_anndata_filename)
        larger_anndata_object = anndata.read_h5ad(self.large_anndata_filename)
        load_dataset(self.corpus_path, small_anndata_object, "dataset_0")
        with tiledb.open(f"{self.corpus_path}/{VAR_ARRAY_NAME}", "r") as var:
            var_df = var.df[:]
            total_stored_genes = len(set(var_df["gene_ontology_term_id"].to_numpy(dtype=str)))
        small_gene_count = 3
        self.assertEqual(small_gene_count, total_stored_genes)

        load_dataset(self.corpus_path, larger_anndata_object, "dataset_1")

        large_gene_count = 1000  # overlapping genes should only be counted once
        with tiledb.open(f"{self.corpus_path}/var", "r") as var:
            var_df = var.df[:]
            total_stored_genes = len(set(var_df["gene_ontology_term_id"].to_numpy(dtype=str)))
        self.assertEqual(large_gene_count, total_stored_genes)

    def test_mapping_between_local_file_and_global_tdb_is_valid_and_consistent_as_datasets_are_added(self):
        """
        DO NOT DELETE THIS TEST
        """
        # load_h5ad(self.path_to_dataset_0, self.corpus_name, False)
        pass

    def test_axes_labels_updated_for_new_genes(self):
        pass

    def test_expression_matrix_matches_global_index(self):
        pass

    def test_raw_expression_matrix_normalized_by_rankit(self):
        pass

    @unittest.skip("removed corpus fixture")
    @patch("backend.wmg.data.cube_pipeline.extract.copy_datasets_to_instance")
    @patch("backend.wmg.data.cube_pipeline.extract.get_dataset_s3_uris")
    def test_corpus_creation_works_as_expected(self, mock_get_uris, mock_copy):
        load_data_and_create_cube(self.path_to_datasets, self.corpus_name, self.tmp_dir)

        # check obs
        with tiledb.open(f"{self.corpus_path}/obs", "r") as obs:
            actual_obs_df = obs.df[:]
        with tiledb.open(self.fixture_file_path("fixtures/small-corpus/obs"), "r") as obs:
            expected_obs_df = obs.df[:]

        self.assertTrue(
            expected_obs_df.development_stage_ontology_term_id.equals(actual_obs_df.development_stage_ontology_term_id)
        )
        self.assertTrue(expected_obs_df.obs_idx.equals(actual_obs_df.obs_idx))
        self.assertTrue(expected_obs_df.dataset_id.equals(actual_obs_df.dataset_id))
        self.assertTrue(expected_obs_df.dataset_local_cell_id.equals(actual_obs_df.dataset_local_cell_id))

        # check vars
        with tiledb.open(f"{self.corpus_path}/var", "r") as var:
            actual_var_df = var.df[:]
        with tiledb.open(self.fixture_file_path("fixtures/small-corpus/var"), "r") as var:
            expected_var_df = var.df[:]

        self.assertTrue(expected_var_df.equals(actual_var_df))

        # check expression matrix
        with tiledb.open(f"{self.corpus_path}/X", "r") as x:
            actual_x_df = x.df[:]
        with tiledb.open(self.fixture_file_path("fixtures/small-corpus/X"), "r") as x:
            expected_x_df = x.df[:]

        self.assertTrue(expected_x_df.equals(actual_x_df))

    def test__filter_out_rankits_with_low_expression_counts__boundaries(self):
        row = [0, 1, 2]
        col = [0, 1, 2]
        rankits = [0.5, 0.7, 0.9]  # 0.5 and 0.7 should be filtered
        raw_counts = [
            RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD - 1,  # should be filtered
            RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD,  # should be filtered
            RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD + 1,  # should not be filtered
        ]
        rankit_csr_matrix = csr_matrix((rankits, (row, col)))
        raw_counts_coo_matrix = coo_matrix((raw_counts, (row, col)))

        rankits_filtered = filter_out_rankits_with_low_expression_counts(rankit_csr_matrix, raw_counts_coo_matrix)

        self.assertEqual(0.9, sum(rankits_filtered.data))

    def test__filter_out_rankits_with_low_expression_counts__majority_filtered(self):
        row = [0, 1]
        col = [0, 1]
        rankits = [0.7, 0.9]  # 0.5 and 0.7 should be filtered
        raw_counts = [
            RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD - 1,  # should be filtered
            RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD + 1,  # should not be filtered
        ]
        rankit_csr_matrix = csr_matrix((rankits, (row, col)))
        raw_counts_coo_matrix = coo_matrix((raw_counts, (row, col)))

        rankits_filtered = filter_out_rankits_with_low_expression_counts(
            rankit_csr_matrix, raw_counts_coo_matrix, expect_majority_filtered=True
        )

        self.assertEqual(0.9, sum(rankits_filtered.data))

    def test__filter_out_rankits_with_low_expression_counts__minority_filtered(self):
        row = [0, 1, 2, 3]
        col = [0, 1, 2, 3]
        rankits = [0.3, 0.5, 0.7, 0.9]  # 0.5 and 0.7 should be filtered
        raw_counts = [
            RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD - 1,  # should be filtered
            RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD + 1,  # should not be filtered
            RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD + 1,  # should not be filtered
            RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD + 1,  # should not be filtered
        ]
        rankit_csr_matrix = csr_matrix((rankits, (row, col)))
        raw_counts_coo_matrix = coo_matrix((raw_counts, (row, col)))

        rankits_filtered = filter_out_rankits_with_low_expression_counts(
            rankit_csr_matrix, raw_counts_coo_matrix, expect_majority_filtered=False
        )

        self.assertEqual(0.5 + 0.7 + 0.9, sum(rankits_filtered.data))

    def test__rankit_scores_ties_the_same(self):
        """
        when the expression values are the same for all genes in the cell the resulting rankit values should all be 3
        """
        expected_rankit_scores = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
        counts = np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]])
        raw_expression_csr_matrix = sparse.csr_matrix(counts)
        rankit_scores = rankit(raw_expression_csr_matrix)
        self.assertEqual(expected_rankit_scores, rankit_scores.data.tolist())

    def test__rankit_value_spread(self):
        """
        Regardless of the expression values, rankit scores should be spread around 3.0 with higher expression levels
        (compared to other gene expression values for the cell) given a larger score and lower expression values given
        a lower score.
        """
        expected_rankit_scores = [2.0325, 3.0, 3.9674, 2.0325, 3.0, 3.9674, 2.0325, 3.0, 3.9674]
        counts = np.array([[1.0, 2.0, 3.0], [1.0, 5.0, 25.0], [1.0, 500.0, 1000.0]])
        raw_expression_csr_matrix = sparse.csr_matrix(counts)
        rankit_scores = rankit(raw_expression_csr_matrix)
        for x in range(len(expected_rankit_scores)):
            expected = expected_rankit_scores[x]
            actual = rankit_scores.data[x]
            self.assertAlmostEqual(expected, actual, 2)

    def test__rankit_handles_zero_values_correctly(self):
        """
        Theoretically there shouldn't be any zero expression values in a sparse matrix
        so they shouldn't receive a rankit score
        """
        expected_rankit_scores = [3.0, 3.0, 3.0, 3.0, 2.325, 3.674]
        counts = np.array([[1.0, 0.0, 1.0], [2.0, 0.0, 2.0], [1.0, 0.0, 3.0]])
        raw_expression_csr_matrix = sparse.csr_matrix(counts)
        rankit_scores = rankit(raw_expression_csr_matrix)
        for x in range(len(expected_rankit_scores)):
            expected = expected_rankit_scores[x]
            actual = rankit_scores.data[x]
            self.assertAlmostEqual(expected, actual, 2)

    @patch("backend.corpus_asset_pipelines.integrated_corpus.load.transform_dataset_raw_counts_to_rankit")
    def test__filter_out_cells_with_incorrect_assays(self, mock_rankit):
        # Create dataset with mixture of included and not included assays
        CELL_COUNT = 5
        test_anndata_object = create_anndata_test_object(num_genes=3, num_cells=CELL_COUNT)
        test_anndata_object.obs["assay_ontology_term_id"].cat.add_categories(["NOT_INCLUDED"], inplace=True)
        test_anndata_object.obs["assay_ontology_term_id"][0] = "NOT_INCLUDED"

        # pre_concatenation filters remove irrelevant data
        updated_test_anndata_object = apply_pre_concatenation_filters(test_anndata_object, min_genes=0)
        load_dataset(self.corpus_path, updated_test_anndata_object, "dataset_0")
        obs = tiledb.open(f"{self.corpus_path}/{OBS_ARRAY_NAME}", "r")
        corpus_cell_count = obs.df[:].shape[0]

        # check the cell count is one less than the starting count
        # because we replaced the assay type for one cell in the original anndata object
        self.assertEqual(corpus_cell_count, CELL_COUNT - 1)

    def test_dataset_validation_checks_correct_expression_matrix(self):
        with self.subTest("Test dataset with sparse X and raw is valid"):
            valid_anndata_object = create_anndata_test_object(num_genes=3, num_cells=5)
            valid_anndata_object.raw = valid_anndata_object
            is_dataset_valid = validate_dataset_properties(valid_anndata_object)
            self.assertTrue(is_dataset_valid)

        with self.subTest("Test dataset with dense X and sparse raw is valid"):
            test_anndata_object = create_anndata_test_object(num_genes=3, num_cells=5)
            dense_matrix = test_anndata_object.X.todense()
            test_anndata_object.raw = test_anndata_object
            test_anndata_object.X = dense_matrix
            is_dataset_valid = validate_dataset_properties(test_anndata_object)
            self.assertTrue(is_dataset_valid)

        with self.subTest("Test dataset with sparse X and dense raw is invalid"):
            test_anndata_object = create_anndata_test_object(num_genes=3, num_cells=5)
            hold_x = test_anndata_object.X
            dense_matrix = test_anndata_object.X.todense()
            test_anndata_object.X = dense_matrix
            test_anndata_object.raw = test_anndata_object
            test_anndata_object.X = hold_x
            is_dataset_valid = validate_dataset_properties(test_anndata_object)
            self.assertFalse(is_dataset_valid)

        with self.subTest("Test dataset with dense X and raw is invalid"):
            invalid_anndata = create_anndata_test_object(num_genes=3, num_cells=5)
            dense_matrix = invalid_anndata.X.todense()
            invalid_anndata.X = dense_matrix
            invalid_anndata.raw = invalid_anndata
            is_dataset_valid = validate_dataset_properties(invalid_anndata)
            self.assertFalse(is_dataset_valid)
