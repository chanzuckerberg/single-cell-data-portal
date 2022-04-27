import os
import pathlib
import shutil
import tempfile
import unittest

import tiledb
from unittest.mock import patch

from scipy.sparse import coo_matrix, csr_matrix

from backend.atlas_asset_pipelines.integrated_corpus.job import build_integrated_corpus, process_h5ad_for_corpus
from backend.atlas_asset_pipelines.integrated_corpus.load import load_h5ad, update_corpus_var, update_corpus_obs
from backend.atlas_asset_pipelines.integrated_corpus.transform import filter_out_rankits_with_low_expression_counts
from backend.wmg.data.constants import RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD
from backend.wmg.data.schemas.corpus_schema import create_tdb, VAR_ARRAY_NAME, OBS_ARRAY_NAME
from tests.unit.backend.wmg.fixtures.test_anndata_object import create_anndata_test_object


class TestCorpusIntegrationETL(unittest.TestCase):
    @staticmethod
    def fixture_file_path(relative_filename):
        parent_dir_of_curr_file = "/".join(os.path.dirname(__file__).split("/")[:-1])
        return os.path.abspath(os.path.join(parent_dir_of_curr_file, relative_filename))

    @classmethod
    def setUpClass(cls) -> None:
        super().setUp(cls)
        cls.tmp_dir = tempfile.mkdtemp()
        cls.BASIC_CELL_COUNT = 5
        cls.BASIC_GENE_COUNT = 3
        cls.LARGER_CELL_COUNT = 5000
        cls.LARGER_GENE_COUNT = 1000
        cls.basic_test_anndata_object = create_anndata_test_object(num_cells=cls.BASIC_CELL_COUNT,
                                                                   num_genes=cls.BASIC_GENE_COUNT)
        cls.larger_test_anndata_object = create_anndata_test_object(num_cells=cls.LARGER_CELL_COUNT,
                                                                    num_genes=cls.LARGER_GENE_COUNT)
        os.mkdir(f"{cls.tmp_dir}/datasets")
        os.mkdir(f"{cls.tmp_dir}/datasets/basic_test_dataset")
        os.mkdir(f"{cls.tmp_dir}/datasets/larger_test_dataset")

        cls.small_anndata_filename = pathlib.Path(cls.tmp_dir, "datasets/basic_test_dataset/local.h5ad")
        cls.large_anndata_filename = pathlib.Path(cls.tmp_dir, "datasets/larger_test_dataset/local.h5ad")

        cls.small_anndata_filename.touch()
        cls.large_anndata_filename.touch()

        cls.basic_test_anndata_object.write(cls.small_anndata_filename, compression="gzip")
        cls.larger_test_anndata_object.write(cls.large_anndata_filename, compression="gzip")

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
            create_tdb(self.tmp_dir, self.corpus_name)
            # self.pbmc3k_anndata_object = get_test_anndata_dataset()

    def tearDown(self) -> None:
        super().tearDown()
        shutil.rmtree(self.corpus_path)

    @patch("backend.atlas_asset_pipelines.integrated_corpus.job.load_h5ad")
    @patch("backend.atlas_asset_pipelines.integrated_corpus.job.validate_dataset_properties")
    def test__datasets_are_validated_before_being_added_to_corpus(self, mock_validation, mock_load):
        mock_validation.return_value = False
        build_integrated_corpus(self.path_to_datasets, self.corpus_path)
        self.assertEqual(mock_load.call_count, 0)

    @patch("backend.atlas_asset_pipelines.integrated_corpus.job.apply_pre_concatenation_filters")
    def test_global_var_array_updated_when_dataset_contains_new_genes(self, mock_filter):
        process_h5ad_for_corpus(self.small_anndata_filename, self.corpus_path, False)
        with tiledb.open(f"{self.corpus_path}/{VAR_ARRAY_NAME}", "r") as var:
            var_df = var.df[:]
            total_stored_genes = len(set(var_df["gene_ontology_term_id"].to_numpy(dtype=str)))
        small_gene_count = 3
        self.assertEqual(small_gene_count, total_stored_genes)

        process_h5ad_for_corpus(self.large_anndata_filename, self.corpus_path, False)

        large_gene_count = 1000  # overlapping genes should only be counted once
        with tiledb.open(f"{self.corpus_path}/{VAR_ARRAY_NAME}", "r") as var:
            var_df = var.df[:]
            total_stored_genes = len(set(var_df["gene_ontology_term_id"].to_numpy(dtype=str)))
        self.assertEqual(large_gene_count, total_stored_genes)

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

    @patch("backend.atlas_asset_pipelines.integrated_corpus.transform.transform_dataset_raw_counts_to_rankit")
    def test__filter_out_cells_with_incorrect_assays(self, mock_rankit):
        # Create dataset with mixture of included and not included assays
        CELL_COUNT = 5
        test_anndata_object = create_anndata_test_object(num_cells=CELL_COUNT, num_genes=3)
        test_anndata_object.obs["assay_ontology_term_id"].cat.add_categories(["NOT_INCLUDED"], inplace=True)
        test_anndata_object.obs["assay_ontology_term_id"][0] = "NOT_INCLUDED"

        os.mkdir(f"{self.tmp_dir}/filter_assays_test_dataset")
        anndata_filename = pathlib.Path(self.tmp_dir, "filter_assays_test_dataset/local.h5ad")
        anndata_filename.touch()
        test_anndata_object.write(anndata_filename, compression="gzip")

        load_h5ad(self.corpus_path, test_anndata_object, dataset_id='lll')
        obs = tiledb.open(f"{self.corpus_path}/obs", "r")
        corpus_cell_count = obs.df[:].shape[0]

        # check the cell count is one less than the starting count
        # because we replaced the assay type for one cell in the original anndata object
        self.assertEqual(corpus_cell_count, CELL_COUNT - 1)

    def test__load_h5ad(self):
        CELL_COUNT = 5
        GENE_COUNT = 3
        test_anndata_object = create_anndata_test_object(num_cells=CELL_COUNT, num_genes=GENE_COUNT)

        os.mkdir(f"{self.tmp_dir}/load_test")
        anndata_filename = pathlib.Path(self.tmp_dir, "load_test/local.h5ad")
        anndata_filename.touch()
        test_anndata_object.write(anndata_filename, compression="gzip")

        load_h5ad(self.corpus_path, test_anndata_object, dataset_id='lll')
        obs = tiledb.open(f"{self.corpus_path}/{OBS_ARRAY_NAME}")
        var = tiledb.open(f"{self.corpus_path}/{VAR_ARRAY_NAME}")

        stored_genes = var.df[:].gene_ontology_term_id
        stored_cells = obs.df[:].dataset_local_cell_id
        self.assertEqual(len(stored_cells), CELL_COUNT)
        self.assertEqual(len(stored_genes), GENE_COUNT)

    def test__update_corpus_var(self):

        initial_corpus_var = update_corpus_var(self.corpus_path, self.basic_test_anndata_object)

        updated_corpus_var = update_corpus_var(self.corpus_path, self.larger_test_anndata_object)

        # check that genes are added
        self.assertEqual(len(initial_corpus_var.index), self.BASIC_GENE_COUNT)

        # check that matching gene ids are only counted once
        self.assertEqual(len(updated_corpus_var.index), self.LARGER_GENE_COUNT)

        # corpus_var indexes on the gene ontology id
        self.assertEqual(updated_corpus_var.index.name, 'gene_ontology_term_id')

    def test__update_corpus_obs_returns_position_of_first_cell(self):
        initial_corpus_obs_starting_idx = update_corpus_obs(self.corpus_path, self.basic_test_anndata_object,
                                                            'basic_test_dataset')
        updated_corpus_obs_starting_idx = update_corpus_obs(self.corpus_path, self.larger_test_anndata_object,
                                                            'larger_test_dataset')
        self.assertEqual(initial_corpus_obs_starting_idx, 0)
        self.assertEqual(updated_corpus_obs_starting_idx, self.BASIC_CELL_COUNT)

    def test__update_corpus_axis_extends_obs_df_by_anndata_obs(self):
        update_corpus_obs(self.corpus_path, self.basic_test_anndata_object,
                          'basic_test_dataset')

        with tiledb.open(f"{self.corpus_path}/{OBS_ARRAY_NAME}", "r") as obs:
            obs_df = obs.df[:]
            cell_id_mapping = obs_df.dataset_local_cell_id

            # map from the corpus index back to the anndata_object index via the cell_id_mapping to ensure data is concatenated correctly
            for i in range(len(cell_id_mapping)):
                self.assertEqual(obs_df.assay_ontology_term_id[i],
                                 self.basic_test_anndata_object.obs[:].assay_ontology_term_id[cell_id_mapping[i]])

        updated_corpus_obs_starting_idx = update_corpus_obs(self.corpus_path, self.larger_test_anndata_object,
                                                            'larger_test_dataset')
        with tiledb.open(f"{self.corpus_path}/{OBS_ARRAY_NAME}", "r") as obs:
            # slice out the data representing the larger_anndata_test_object obs
            larger_test_anndata_object_obs = obs.df[:][updated_corpus_obs_starting_idx:]
            cell_id_mapping = larger_test_anndata_object_obs.dataset_local_cell_id
            for i in range(updated_corpus_obs_starting_idx, len(cell_id_mapping)):
                self.assertEqual(larger_test_anndata_object_obs.assay_ontology_term_id[i],
                                 self.larger_test_anndata_object.obs[:].assay_ontology_term_id[cell_id_mapping[i]])

    @unittest.skip("TO DO IMPLEMENT")
    def test_var_labels_contain_correct_dimensions(self):
        raise NotImplementedError

    @unittest.skip("TO DO IMPLEMENT")
    def test_label_type(self):
        raise NotImplementedError
