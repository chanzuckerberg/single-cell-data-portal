import json
import unittest

from backend.cellguide.pipeline.computational_marker_genes import get_computational_marker_genes
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils.dict_compare import compare_dicts
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)
from tests.unit.cellguide_pipeline.constants import (
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    COMPUTATIONAL_MARKER_GENES_FIXTURE_FILENAME,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class MarkerGeneCalculatorTests(unittest.TestCase):
    def test__marker_gene_calculation(self):
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{COMPUTATIONAL_MARKER_GENES_FIXTURE_FILENAME}", "r") as f:
            expected__computational_marker_genes = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            cell_counts_df = snapshot.cell_counts_cube.df[:]
            tree_builder = OntologyTreeBuilder(cell_counts_df)
            computational_marker_genes = convert_dataclass_to_dict_and_strip_nones(
                get_computational_marker_genes(
                    snapshot=snapshot,
                    ontology_tree=tree_builder,
                )
            )
            self.assertTrue(compare_dicts(computational_marker_genes, expected__computational_marker_genes))
