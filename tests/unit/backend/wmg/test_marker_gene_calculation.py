import unittest
import json
import pytest
from tests.unit.backend.wmg.fixtures.test_snapshot import load_test_fmg_snapshot
from backend.wmg.pipeline.summary_cubes.calculate_markers import _query_tiledb, get_markers
from backend.wmg.data.query import retrieve_top_n_markers

TEST_SNAPSHOT = "test-fmg-snapshot"

TARGET_FILTERS = {
    "tissue_ontology_term_ids": ["UBERON:0002048"],
    "cell_type_ontology_term_ids": ["CL:0000786"],
    "organism_ontology_term_id": "NCBITaxon:9606",
    "disease_ontology_term_ids": ["PATO:0000461"],
}
CONTEXT_FILTERS = {
    "tissue_ontology_term_ids": ["UBERON:0002048"],
    "organism_ontology_term_id": "NCBITaxon:9606",
    "disease_ontology_term_ids": ["PATO:0000461"],
}

# tests key functions in calculate_markers.py in order from most to least nested


class MarkerGeneCalculationTest(unittest.TestCase):
    def test__query_tiledb(self):
        with load_test_fmg_snapshot(TEST_SNAPSHOT) as snapshot:
            agg, t_n_cells_sum, _ = _query_tiledb(TARGET_FILTERS, corpus=snapshot)
            test_sum = list(agg.sum(0))
            # check that returned dataframe is correct
            expected_sum = [21988.90625, 66193.2578125, 8782.0, 8666.0]
            for i in range(len(test_sum)):
                assert abs(test_sum[i] - expected_sum[i]) < 0.05

            # check that returned population sizes are correct
            assert sum(list(t_n_cells_sum.values())[0]) == 70234

    def test__get_markers_ttest(self):
        with load_test_fmg_snapshot(TEST_SNAPSHOT) as snapshot:
            result = get_markers(TARGET_FILTERS, CONTEXT_FILTERS, corpus=snapshot, test="ttest", percentile=0.05)
            result = json.loads(
                json.dumps(result).replace("p_value_ttest", "p_value").replace("effect_size_ttest", "effect_size")
            )
            result = [
                {"gene_ontology_term_id": k, "p_value": v["p_value"], "effect_size": v["effect_size"]}
                for k, v in result.items()
            ]
            tissue = TARGET_FILTERS["tissue_ontology_term_ids"]
            celltype = TARGET_FILTERS["cell_type_ontology_term_ids"]
            organism = TARGET_FILTERS["organism_ontology_term_id"]

            expected = snapshot.marker_genes_cube.df[(tissue, organism, celltype)]
            expected = retrieve_top_n_markers(expected, "ttest", 10)
            for i, elem in enumerate(result):
                assert pytest.approx(elem) == expected[i]

    def test__get_markers_binomtest(self):
        with load_test_fmg_snapshot(TEST_SNAPSHOT) as snapshot:
            result = get_markers(TARGET_FILTERS, CONTEXT_FILTERS, corpus=snapshot, test="binomtest", percentile=0.3)
            result = json.loads(
                json.dumps(result)
                .replace("p_value_binomtest", "p_value")
                .replace("effect_size_binomtest", "effect_size")
            )
            result = [
                {"gene_ontology_term_id": k, "p_value": v["p_value"], "effect_size": v["effect_size"]}
                for k, v in result.items()
            ]
            tissue = TARGET_FILTERS["tissue_ontology_term_ids"]
            celltype = TARGET_FILTERS["cell_type_ontology_term_ids"]
            organism = TARGET_FILTERS["organism_ontology_term_id"]

            expected = snapshot.marker_genes_cube.df[(tissue, organism, celltype)]
            expected = retrieve_top_n_markers(expected, "binomtest", 10)
            for i, elem in enumerate(result):
                assert pytest.approx(elem) == expected[i]
