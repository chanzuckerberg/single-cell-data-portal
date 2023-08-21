import json
import unittest

from backend.cellguide.pipeline.computational_marker_genes import get_marker_genes_per_and_across_tissues
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils.dict_compare import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class MarkerGeneCalculatorTests(unittest.TestCase):
    def test__marker_gene_calculation(self):
        with open("tests/unit/cellguide_pipeline/fixtures/computational_marker_genes.json", "r") as f:
            expected__computational_marker_genes = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            cell_counts_df = snapshot.cell_counts_cube.df[:]
            tree_builder = OntologyTreeBuilder(cell_counts_df)
            computational_marker_genes = convert_dataclass_to_dict_and_strip_nones(
                get_marker_genes_per_and_across_tissues(
                    snapshot=snapshot,
                    all_cell_types_in_corpus=tree_builder.all_cell_type_ids_in_corpus,
                )
            )

            self.assertTrue(compare_dicts(computational_marker_genes, expected__computational_marker_genes))
