import unittest
from unittest import mock

from backend.wmg.data.cube_pipeline import load
from tests.unit.backend.corpora.fixtures.environment_setup import fixture_file_path


class TestCorpusLoad(unittest.TestCase):
    def setUp(self) -> None:
        self.path_to_dataset_0 = fixture_file_path('small_datasets/be573215-4116-40a1-9c0f-1b21ba53482b/local.h5ad')
        self.path_to_dataset_1 = fixture_file_path('small_datasets/d50b8959-6ce9-4a9b-b804-99892c93b183/local.h5ad')
        self.path_to_datasets = fixture_file_path('small_datasets')

    @mock.patch('backend.wmg.load_corpus.load_h5ad')
    def test__load_loads_all_datasets_in_directory(self, mock_load_h5ad):
        load(self.path_to_datasets, 'test-group')
        import pdb
        pdb.set_trace()
        self.assertCalled

    @mock.patch('backend.wmg.load_corpus.load_h5ad')
    def test_loading_a_dataset_already_in_the_corpus_is_a_noop(self):
        pass

    def test_sparse_datasets_are_not_added_to_corpus(self):
        pass

    def test_datasets_with_the_wrong_schema_are_not_added_to_corpus(self):
        pass

    def test_global_var_array_updated_when_dataset_contains_new_genes(self):
        pass

    def test_mapping_between_local_file_and_global_tdb_is_valid_and_consistent_as_datasets_are_added(self):
        """
        DO NOT DELETE THIS TEST
        """
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

    def test_raw_expression_matrix_normalized_by_rankit():
        pass
