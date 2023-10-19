import tempfile
import unittest

import tiledb

from backend.wmg.data.schemas.cube_schema import cell_counts_schema, expression_summary_schema
from backend.wmg.pipeline.utils import create_empty_cube_if_needed


class TestSummaryCubeSchema(unittest.TestCase):
    def test_summary_cube_schema_contains_correct_dimensions(self):
        with tempfile.TemporaryDirectory() as summary_cube_dir:
            create_empty_cube_if_needed(uri=f"{summary_cube_dir}/cube", schema=expression_summary_schema)
            summary_cube_schema = tiledb.ArraySchema.load(f"{summary_cube_dir}/cube")
            self.assertEqual(summary_cube_schema.ndim, 3)
            self.assertTrue(summary_cube_schema.domain.has_dim("gene_ontology_term_id"))
            self.assertTrue(summary_cube_schema.domain.has_dim("tissue_ontology_term_id"))
            self.assertTrue(summary_cube_schema.domain.has_dim("organism_ontology_term_id"))

    def test_summary_cube_schema_contains_correct_attrs(self):
        with tempfile.TemporaryDirectory() as summary_cube_dir:
            create_empty_cube_if_needed(uri=f"{summary_cube_dir}/cube", schema=expression_summary_schema)
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


class TestCellCountCubeSchema(unittest.TestCase):
    def test_cell_count_cube_schema_contains_correct_dimensions(self):
        with tempfile.TemporaryDirectory() as cell_count_cube_dir:
            create_empty_cube_if_needed(uri=f"{cell_count_cube_dir}/cube", schema=cell_counts_schema)
            cell_count_cube_schema = tiledb.ArraySchema.load(f"{cell_count_cube_dir}/cube")
            self.assertEqual(cell_count_cube_schema.ndim, 2)
            self.assertTrue(cell_count_cube_schema.domain.has_dim("tissue_ontology_term_id"))
            self.assertTrue(cell_count_cube_schema.domain.has_dim("organism_ontology_term_id"))

    def test_cell_count_cube_schema_contains_correct_attributes(self):
        with tempfile.TemporaryDirectory() as cell_count_cube_dir:
            create_empty_cube_if_needed(uri=f"{cell_count_cube_dir}/cube", schema=cell_counts_schema)
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
