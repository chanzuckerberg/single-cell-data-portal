import os
import pathlib
import shutil
import tempfile
import unittest
from unittest.mock import Mock, patch

import anndata
import numpy as np
import tiledb
from scipy import sparse
from scipy.sparse import coo_matrix, csr_matrix

from backend.wmg.data.constants import RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD
from backend.wmg.data.rankit import rankit
from backend.wmg.data.schemas.corpus_schema import OBS_ARRAY_NAME, VAR_ARRAY_NAME, create_tdb_integrated_corpus
from backend.wmg.pipeline.cube_pipeline import load_data_and_create_cube
from backend.wmg.pipeline.integrated_corpus.job import build_integrated_corpus
from backend.wmg.pipeline.integrated_corpus.load import load_dataset
from backend.wmg.pipeline.integrated_corpus.transform import (
    apply_pre_concatenation_filters,
    filter_out_rankits_with_low_expression_counts,
)
from backend.wmg.pipeline.integrated_corpus.validate import validate_dataset_properties
from tests.unit.backend.wmg.fixtures.test_anndata_object import create_anndata_test_fixture, create_anndata_test_object


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
        cls.path_to_datasets = pathlib.Path(cls.tmp_dir, "datasets")
        os.mkdir(cls.path_to_datasets)
        cls.small_anndata_filename = create_anndata_test_fixture(cls.path_to_datasets, "basic_test_dataset", 3, 5)
        cls.large_anndata_filename = create_anndata_test_fixture(cls.path_to_datasets, "large_test_dataset", 1000, 5000)

    @classmethod
    def tearDownClass(cls) -> None:
        super().tearDownClass()
        shutil.rmtree(cls.tmp_dir)

    def setUp(self) -> None:
        super().setUp()
        self.path_to_datasets = pathlib.Path(self.tmp_dir, "datasets")
        self.corpus_path = f"{self.tmp_dir}/test-corpus"
        if not tiledb.VFS().is_dir(self.corpus_path):
            create_tdb_integrated_corpus(self.corpus_path)
            # self.pbmc3k_anndata_object = get_test_anndata_dataset()

    def tearDown(self) -> None:
        super().tearDown()
        shutil.rmtree(self.corpus_path)

    @patch("backend.wmg.pipeline.integrated_corpus.job.tiledb.consolidate")
    @patch("backend.wmg.pipeline.integrated_corpus.job.tiledb.vacuum")
    @patch("backend.wmg.pipeline.integrated_corpus.job.process_h5ad_for_corpus")
    def test__load_loads_all_datasets_in_directory(self, mock_process_h5ad, mock_vacuum, mock_consolidate):
        mock_process_h5ad.return_value = None, None
        build_integrated_corpus(self.path_to_datasets, self.corpus_path)
        self.assertEqual(mock_process_h5ad.call_count, 2)
        self.assertEqual(mock_vacuum.call_count, 3)
        self.assertEqual(mock_consolidate.call_count, 3)

    @patch("backend.wmg.pipeline.integrated_corpus.load.update_corpus_var")
    @patch("backend.wmg.pipeline.integrated_corpus.job.validate_dataset_properties")
    def test_invalid_datasets_are_not_added_to_corpus(self, mock_validation, mock_global_var):
        mock_validation.return_value = False
        build_integrated_corpus(self.path_to_datasets, self.corpus_path)
        self.assertEqual(mock_global_var.call_count, 0)

    @patch("backend.wmg.pipeline.integrated_corpus.load.transform_dataset_raw_counts_to_rankit")
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

    def test_axes_labels_updated_for_new_genes(self):
        pass

    def test_expression_matrix_matches_global_index(self):
        pass

    def test_raw_expression_matrix_normalized_by_rankit(self):
        pass

    @patch("backend.wmg.pipeline.integrated_corpus.transform.GENE_EXPRESSION_COUNT_MIN_THRESHOLD", 1)
    @patch("backend.wmg.pipeline.integrated_corpus.job.tiledb.consolidate", new=Mock())  # Slow
    @patch("backend.wmg.pipeline.integrated_corpus.job.tiledb.vacuum", new=Mock())  # Slow
    @patch("backend.wmg.pipeline.cube_pipeline.upload_artifacts_to_s3", new=Mock())  # Don't upload the cube.
    @patch("backend.wmg.pipeline.integrated_corpus.job.extract.get_dataset_asset_urls", new=Mock(return_value={}))
    @unittest.skip("Temporarily skipping until we can resolve cell counts cube bug")
    def test_snapshot_creation_works_as_expected(self):
        generate_cells = 5000
        expected_datasets = 2
        expected_genes = 1000
        expected_cell_count = 0

        with tempfile.TemporaryDirectory() as temp_dir:
            path_to_datasets = f"{temp_dir}/datasets"
            os.mkdir(path_to_datasets)
            for i in range(expected_datasets):
                create_anndata_test_fixture(path_to_datasets, f"dataset_{i}", expected_genes, generate_cells)
                expected_cell_count = expected_cell_count + generate_cells

            # Run
            snapshot_path, stats = load_data_and_create_cube(path_to_datasets, self.tmp_dir, validate_cube=False)

            # Verify
            self.assertEqual(stats["cell_count"], expected_cell_count)
            self.assertEqual(stats["gene_count"], expected_genes)
            self.assertEqual(stats["dataset_count"], expected_datasets)

            # check obs matrix
            with tiledb.open(f"{snapshot_path}/obs", "r") as obs:
                obs_df = obs.df[:]

                # total number of rows should be the number of cells generated per dataset * number of datasets
                assert len(obs_df) == generate_cells * expected_datasets

                # two expected datasets: dataset_0 and dataset_1
                assert (obs_df.obs_idx >= 0).all()
                assert sorted(obs_df.dataset_id.unique()) == ["dataset_0", "dataset_1"]

                # cell ids should be between Cell_0 and Cell_4999
                assert (obs_df.dataset_local_cell_id.str.startswith("Cell_")).all()
                assert obs_df.dataset_local_cell_id.nunique() == generate_cells

                # these are fixed values indicating lung tissue
                assert (obs_df.tissue_original_ontology_term_id == "UBERON:0000101").all()
                assert (obs_df.tissue_ontology_term_id == "UBERON:0002048").all()

            # check var matrix
            with tiledb.open(f"{self.corpus_path}/var", "r") as var:
                var_df = var.df[:]
                assert len(var_df) == 0
                expected_columns = ["feature_name", "feature_reference", "gene_ontology_term_id", "var_idx"]
                assert sorted(var_df.columns) == expected_columns

            # check integrated expression matrix
            with tiledb.open(f"{self.corpus_path}/integrated", "r") as ie:
                ie_df = ie.df[:]
                assert len(var_df) == 0
                expected_columns = ["obs_idx", "rankit", "var_idx"]
                assert sorted(ie_df.columns) == expected_columns

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

    @patch("backend.wmg.pipeline.integrated_corpus.load.transform_dataset_raw_counts_to_rankit")
    def test__filter_out_cells_with_incorrect_assays(self, mock_rankit):
        # Create dataset with mixture of included and not included assays
        CELL_COUNT = 5
        test_anndata_object = create_anndata_test_object(num_genes=3, num_cells=CELL_COUNT)
        test_anndata_object.obs["assay_ontology_term_id"] = test_anndata_object.obs[
            "assay_ontology_term_id"
        ].cat.add_categories(["NOT_INCLUDED"])
        test_anndata_object.obs["assay_ontology_term_id"][0] = "NOT_INCLUDED"

        # pre_concatenation filters remove irrelevant data
        updated_test_anndata_object = apply_pre_concatenation_filters(test_anndata_object, min_genes=0)
        load_dataset(self.corpus_path, updated_test_anndata_object, "dataset_0")
        obs_array = tiledb.open(f"{self.corpus_path}/{OBS_ARRAY_NAME}", "r")

        obs_df = obs_array.df[:]
        obs_after_filtering = obs_df[obs_df["filter_cells"] is False]
        corpus_cell_count = obs_after_filtering.shape[0]

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
