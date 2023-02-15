import json
import unittest

import pytest

from backend.wmg.data.query import retrieve_top_n_markers
from backend.wmg.pipeline.summary_cubes.calculate_markers import _prepare_indices_and_metrics, get_markers
from tests.unit.backend.wmg.fixtures.test_snapshot import load_test_fmg_snapshot

TEST_SNAPSHOT = "test-fmg-snapshot"

TARGET_FILTERS = {
    "tissue_ontology_term_ids": ["UBERON:0002048"],
    "cell_type_ontology_term_ids": ["CL:0000786"],
    "organism_ontology_term_id": "NCBITaxon:9606",
}
CONTEXT_FILTERS = {
    "tissue_ontology_term_ids": ["UBERON:0002048"],
    "organism_ontology_term_id": "NCBITaxon:9606",
}

# tests key functions in calculate_markers.py in order from most to least nested


class MarkerGeneCalculationTest(unittest.TestCase):
    def test__query_tiledb(self):
        with load_test_fmg_snapshot(TEST_SNAPSHOT) as snapshot:
            output = _prepare_indices_and_metrics(TARGET_FILTERS, CONTEXT_FILTERS, corpus=snapshot)
            context_agg = output[0]
            target_agg = output[1]
            n_cells_per_gene_target = output[3]
            n_cells_per_gene_context = output[4]

            test_sum_target = list(target_agg.sum(0))
            # check that returned dataframe is correct
            expected_sum_target = [34725.328125, 103803.515625, 13914.0, 13777.0]
            for i in range(len(test_sum_target)):
                assert abs(test_sum_target[i] - expected_sum_target[i]) < 0.05

            test_sum_context = list(context_agg.sum(0))
            # check that returned dataframe is correct
            expected_sum_context = [25811448.0, 54070508.0, 14401558.0, 13712079.0]
            for i in range(len(test_sum_context)):
                assert abs(test_sum_context[i] - expected_sum_context[i]) < 0.05

            # check that returned population sizes are correct
            assert n_cells_per_gene_target.sum() == 302270.0
            assert n_cells_per_gene_context.sum() == 901022918.0

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
