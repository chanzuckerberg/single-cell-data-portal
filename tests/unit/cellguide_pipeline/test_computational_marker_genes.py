import json
import unittest
from unittest.mock import patch

from backend.cellguide.pipeline.computational_marker_genes import get_computational_marker_genes
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import convert_dataclass_to_dict_and_strip_nones
from tests.test_utils import compare_dicts
from tests.test_utils.mocks import mock_bootstrap_rows_percentiles
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)
from tests.unit.cellguide_pipeline.constants import (
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    COMPUTATIONAL_MARKER_GENES_FIXTURE_FILENAME,
    REFORMATTED_COMPUTATIONAL_MARKER_GENES_FIXTURE_FILENAME,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class MarkerGeneCalculatorTests(unittest.TestCase):
    @unittest.skip("This test fails in GHA but passes locally - skipping until it is fixed (issue #6240).")
    def test__marker_gene_calculation(self):
        with open(f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{COMPUTATIONAL_MARKER_GENES_FIXTURE_FILENAME}", "r") as f:
            expected__computational_marker_genes = json.load(f)
        with open(
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{REFORMATTED_COMPUTATIONAL_MARKER_GENES_FIXTURE_FILENAME}", "r"
        ) as f:
            expected__reformatted_marker_genes = json.load(f)
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot, patch(
            "backend.cellguide.pipeline.computational_marker_genes.computational_markers.bootstrap_rows_percentiles",
            new=mock_bootstrap_rows_percentiles,
        ):
            cell_counts_df = snapshot.cell_counts_cube.df[:]
            tree_builder = OntologyTreeBuilder(cell_counts_df)
            computational_marker_genes, reformatted_marker_genes = get_computational_marker_genes(
                snapshot=snapshot,
                ontology_tree=tree_builder,
            )
            computational_marker_genes = convert_dataclass_to_dict_and_strip_nones(computational_marker_genes)
            reformatted_marker_genes = convert_dataclass_to_dict_and_strip_nones(reformatted_marker_genes)
            self.assertTrue(compare_dicts(computational_marker_genes, expected__computational_marker_genes))
            self.assertTrue(compare_dicts(reformatted_marker_genes, expected__reformatted_marker_genes))
