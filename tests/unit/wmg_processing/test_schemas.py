import tempfile
import unittest

import tiledb

from backend.wmg.data.schemas.corpus_schema import create_tdb_integrated_corpus
from backend.wmg.data.schemas.cube_schema import cell_counts_schema, expression_summary_schema
from backend.wmg.data.schemas.expression_summary_fmg_cube_schema import expression_summary_fmg_schema
from backend.wmg.data.utils import create_empty_cube


class TestIntegratedCorpusSchema(unittest.TestCase):
    def test_integrated_corpus_var_schema_contain_correct_dimensions(self):
        with tempfile.TemporaryDirectory() as corpus_dir:
            corpus_path = f"{corpus_dir}/test_corpus"
            create_tdb_integrated_corpus(corpus_path)
            var_schema = tiledb.ArraySchema.load(f"{corpus_path}/var")
            self.assertEqual(var_schema.domain.ndim, 1)
            self.assertTrue(var_schema.domain.has_dim("gene_ontology_term_id"))

    def test_integrated_corpus_var_schema_contain_correct_attributes(self):
        with tempfile.TemporaryDirectory() as corpus_dir:
            corpus_path = f"{corpus_dir}/test_corpus"
            create_tdb_integrated_corpus(corpus_path)
            var_schema = tiledb.ArraySchema.load(f"{corpus_path}/var")
            self.assertEqual(var_schema.nattr, 3)
            self.assertTrue(var_schema.has_attr("var_idx"))
            self.assertTrue(var_schema.has_attr("feature_reference"))
            self.assertTrue(var_schema.has_attr("feature_name"))

    def test_integrated_corpus_obs_schema_contain_correct_dimensions(self):
        with tempfile.TemporaryDirectory() as corpus_dir:
            corpus_path = f"{corpus_dir}/test_corpus"
            create_tdb_integrated_corpus(corpus_path)
            obs_schema = tiledb.ArraySchema.load(f"{corpus_path}/obs")
            self.assertEqual(obs_schema.domain.ndim, 4)
            self.assertTrue(obs_schema.domain.has_dim("dataset_id"))
            self.assertTrue(obs_schema.domain.has_dim("cell_type_ontology_term_id"))
            self.assertTrue(obs_schema.domain.has_dim("tissue_ontology_term_id"))

    def test_integrated_corpus_obs_schema_contain_correct_attributes(self):
        with tempfile.TemporaryDirectory() as corpus_dir:
            corpus_path = f"{corpus_dir}/test_corpus"
            create_tdb_integrated_corpus(corpus_path)
            obs_schema = tiledb.ArraySchema.load(f"{corpus_path}/obs")
            self.assertEqual(obs_schema.nattr, 17)
            self.assertTrue(obs_schema.has_attr("obs_idx"))
            self.assertTrue(obs_schema.has_attr("filter_cells"))
            self.assertTrue(obs_schema.has_attr("cell_type"))
            self.assertTrue(obs_schema.has_attr("assay"))
            self.assertTrue(obs_schema.has_attr("assay_ontology_term_id"))
            self.assertTrue(obs_schema.has_attr("development_stage"))
            self.assertTrue(obs_schema.has_attr("development_stage_ontology_term_id"))
            self.assertTrue(obs_schema.has_attr("disease_ontology_term_id"))
            self.assertTrue(obs_schema.has_attr("tissue"))
            self.assertTrue(obs_schema.has_attr("self_reported_ethnicity"))
            self.assertTrue(obs_schema.has_attr("self_reported_ethnicity_ontology_term_id"))
            self.assertTrue(obs_schema.has_attr("sex"))
            self.assertTrue(obs_schema.has_attr("sex_ontology_term_id"))
            self.assertTrue(obs_schema.has_attr("organism"))
            self.assertTrue(obs_schema.has_attr("organism_ontology_term_id"))
            self.assertTrue(obs_schema.has_attr("dataset_local_cell_id"))


class TestSummaryCubeSchema(unittest.TestCase):
    def test_summary_cube_schema_contains_correct_dimensions(self):
        with tempfile.TemporaryDirectory() as summary_cube_dir:
            create_empty_cube(uri=f"{summary_cube_dir}/cube", schema=expression_summary_schema)
            summary_cube_schema = tiledb.ArraySchema.load(f"{summary_cube_dir}/cube")
            self.assertEqual(summary_cube_schema.ndim, 3)
            self.assertTrue(summary_cube_schema.domain.has_dim("gene_ontology_term_id"))
            self.assertTrue(summary_cube_schema.domain.has_dim("tissue_ontology_term_id"))
            self.assertTrue(summary_cube_schema.domain.has_dim("organism_ontology_term_id"))

    def test_summary_cube_schema_contains_correct_attrs(self):
        with tempfile.TemporaryDirectory() as summary_cube_dir:
            create_empty_cube(uri=f"{summary_cube_dir}/cube", schema=expression_summary_schema)
            summary_cube_schema = tiledb.ArraySchema.load(f"{summary_cube_dir}/cube")
            self.assertEqual(summary_cube_schema.nattr, 11)
            self.assertTrue(summary_cube_schema.has_attr("cell_type_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("tissue_original_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("dataset_id"))
            self.assertTrue(summary_cube_schema.has_attr("assay_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("publication_citation"))
            self.assertTrue(summary_cube_schema.has_attr("development_stage_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("disease_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("self_reported_ethnicity_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("sex_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("nnz"))
            self.assertTrue(summary_cube_schema.has_attr("sum"))


class TestFmgSummaryCubeSchema(unittest.TestCase):
    def test_fmg_summary_cube_schema_contains_correct_dimensions(self):
        with tempfile.TemporaryDirectory() as summary_cube_dir:
            create_empty_cube(uri=f"{summary_cube_dir}/cube", schema=expression_summary_fmg_schema)
            summary_cube_schema = tiledb.ArraySchema.load(f"{summary_cube_dir}/cube")
            self.assertEqual(summary_cube_schema.ndim, 3)
            self.assertTrue(summary_cube_schema.domain.has_dim("cell_type_ontology_term_id"))
            self.assertTrue(summary_cube_schema.domain.has_dim("tissue_ontology_term_id"))
            self.assertTrue(summary_cube_schema.domain.has_dim("organism_ontology_term_id"))

    def test_fmg_summary_cube_schema_contains_correct_attrs(self):
        with tempfile.TemporaryDirectory() as summary_cube_dir:
            create_empty_cube(uri=f"{summary_cube_dir}/cube", schema=expression_summary_fmg_schema)
            summary_cube_schema = tiledb.ArraySchema.load(f"{summary_cube_dir}/cube")
            self.assertEqual(summary_cube_schema.nattr, 12)
            self.assertTrue(summary_cube_schema.has_attr("gene_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("tissue_original_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("dataset_id"))
            self.assertTrue(summary_cube_schema.has_attr("assay_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("development_stage_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("disease_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("self_reported_ethnicity_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("sex_ontology_term_id"))
            self.assertTrue(summary_cube_schema.has_attr("nnz"))
            self.assertTrue(summary_cube_schema.has_attr("sum"))
            self.assertTrue(summary_cube_schema.has_attr("sqsum"))
            self.assertTrue(summary_cube_schema.has_attr("nnz_thr"))


class TestCellCountCubeSchema(unittest.TestCase):
    def test_cell_count_cube_schema_contains_correct_dimensions(self):
        with tempfile.TemporaryDirectory() as cell_count_cube_dir:
            create_empty_cube(uri=f"{cell_count_cube_dir}/cube", schema=cell_counts_schema)
            cell_count_cube_schema = tiledb.ArraySchema.load(f"{cell_count_cube_dir}/cube")
            self.assertEqual(cell_count_cube_schema.ndim, 2)
            self.assertTrue(cell_count_cube_schema.domain.has_dim("tissue_ontology_term_id"))
            self.assertTrue(cell_count_cube_schema.domain.has_dim("organism_ontology_term_id"))

    def test_cell_count_cube_schema_contains_correct_attributes(self):
        with tempfile.TemporaryDirectory() as cell_count_cube_dir:
            create_empty_cube(uri=f"{cell_count_cube_dir}/cube", schema=cell_counts_schema)
            cell_count_cube_schema = tiledb.ArraySchema.load(f"{cell_count_cube_dir}/cube")
            self.assertEqual(cell_count_cube_schema.nattr, 10)
            self.assertTrue(cell_count_cube_schema.has_attr("cell_type_ontology_term_id"))
            self.assertTrue(cell_count_cube_schema.has_attr("tissue_original_ontology_term_id"))
            self.assertTrue(cell_count_cube_schema.has_attr("dataset_id"))
            self.assertTrue(cell_count_cube_schema.has_attr("assay_ontology_term_id"))
            self.assertTrue(cell_count_cube_schema.has_attr("development_stage_ontology_term_id"))
            self.assertTrue(cell_count_cube_schema.has_attr("publication_citation"))
            self.assertTrue(cell_count_cube_schema.has_attr("disease_ontology_term_id"))
            self.assertTrue(cell_count_cube_schema.has_attr("self_reported_ethnicity_ontology_term_id"))
            self.assertTrue(cell_count_cube_schema.has_attr("sex_ontology_term_id"))
            self.assertTrue(cell_count_cube_schema.has_attr("n_cells"))
