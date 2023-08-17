import json
import unittest

import pytest

from backend.wmg.data.query import retrieve_top_n_markers
from backend.wmg.pipeline.summary_cubes.calculate_markers import get_markers
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot

TEST_SNAPSHOT = "realistic-test-snapshot"

TARGET_FILTERS = {
    "tissue_ontology_term_ids": ["UBERON:0002048"],
    "cell_type_ontology_term_ids": ["CL:0000786"],
    "organism_ontology_term_id": "NCBITaxon:9606",
}
CONTEXT_FILTERS = {
    "tissue_ontology_term_ids": ["UBERON:0002048"],
    "organism_ontology_term_id": "NCBITaxon:9606",
}

# Tests key functions in calculate_markers.py in order from most to least nested


class MarkerGeneCalculationTest(unittest.TestCase):
    def test__get_markers_ttest(self):
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            celltype = "CL:0000786"
            tissue = "UBERON:0002048"
            organism = "NCBITaxon:9606"
            result = get_markers(celltype, tissue, organism, corpus=snapshot, percentile=0.05)
            result = json.loads(
                json.dumps(result).replace("p_value_ttest", "p_value").replace("effect_size_ttest", "effect_size")
            )
            result = [
                {"gene_ontology_term_id": k, "p_value": v["p_value"], "effect_size": v["effect_size"]}
                for k, v in result.items()
            ]

            expected = snapshot.marker_genes_cube.df[(tissue, organism, celltype)]
            expected = retrieve_top_n_markers(expected, "ttest", 10)
            for i, elem in enumerate(result):
                assert pytest.approx(elem, rel=0.01) == expected[i]
