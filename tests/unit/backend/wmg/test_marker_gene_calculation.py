import unittest
from tests.unit.backend.wmg.fixtures.test_snapshot import load_test_fmg_snapshot
from backend.wmg.api.calculate_markers import _query_tiledb, get_markers

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
            agg, t_n_cells_sum, _ = _query_tiledb(TARGET_FILTERS, corpus=snapshot)
            test_sum = list(agg.sum(0))
            # check that returned dataframe is correct
            expected_sum = [28538.255859375, 85875.046875, 11312.0, 11185.0]
            for i in range(len(test_sum)):
                assert abs(test_sum[i] - expected_sum[i]) < 0.05

            # check that returned population sizes are correct
            assert sum(list(t_n_cells_sum.values())[0]) == 98974

    def test__get_markers_ttest(self):
        with load_test_fmg_snapshot(TEST_SNAPSHOT) as snapshot:
            result = get_markers(TARGET_FILTERS, CONTEXT_FILTERS, corpus=snapshot, test="ttest")

            expected = {
                "ENSG00000170476": {"p_value_ttest": 2.853570751404757e-277, "effect_size_ttest": 2.8502310053712137},
                "ENSG00000180879": {"p_value_ttest": 3.2048802928401952e-208, "effect_size_ttest": 2.5837366521752467},
                "ENSG00000132465": {"p_value_ttest": 7.223446844565219e-271, "effect_size_ttest": 2.3984059302087104},
                "ENSG00000134285": {"p_value_ttest": 2.4746219866006915e-227, "effect_size_ttest": 1.8600029299922847},
                "ENSG00000166562": {"p_value_ttest": 1.9153063397143112e-215, "effect_size_ttest": 1.8400333201504637},
                "ENSG00000051108": {"p_value_ttest": 2.842480381188694e-211, "effect_size_ttest": 1.6856227945515443},
                "ENSG00000099958": {"p_value_ttest": 3.357622968333614e-205, "effect_size_ttest": 1.6013719201864456},
                "ENSG00000211592": {"p_value_ttest": 7.105551700160005e-93, "effect_size_ttest": 1.4867776488396038},
                "ENSG00000118363": {"p_value_ttest": 1.7808284501343359e-183, "effect_size_ttest": 1.476018358457155},
                "ENSG00000100219": {"p_value_ttest": 1.2807063177018397e-178, "effect_size_ttest": 1.4613032640890151},
            }
            self.assertDictEqual(result, expected)

    def test__get_markers_binomtest(self):
        with load_test_fmg_snapshot(TEST_SNAPSHOT) as snapshot:
            result = get_markers(TARGET_FILTERS, CONTEXT_FILTERS, corpus=snapshot, test="binomtest")
            expected = {
                "ENSG00000099958": {
                    "p_value_binomtest": 1.9613120062451643e-147,
                    "effect_size_binomtest": 7.27287154339715,
                },
                "ENSG00000170476": {
                    "p_value_binomtest": 7.046819867598114e-151,
                    "effect_size_binomtest": 5.982245009212843,
                },
                "ENSG00000134285": {
                    "p_value_binomtest": 8.680479676114522e-137,
                    "effect_size_binomtest": 3.6003438891602255,
                },
                "ENSG00000132465": {
                    "p_value_binomtest": 2.7014577109759077e-141,
                    "effect_size_binomtest": 3.5873641915459094,
                },
                "ENSG00000166562": {
                    "p_value_binomtest": 4.700596957883503e-129,
                    "effect_size_binomtest": 2.3480413711693804,
                },
                "ENSG00000173110": {
                    "p_value_binomtest": 5.911640500830793e-28,
                    "effect_size_binomtest": 1.802670539695546,
                },
                "ENSG00000211592": {
                    "p_value_binomtest": 8.619975795729593e-90,
                    "effect_size_binomtest": 1.563336930583122,
                },
                "ENSG00000171311": {
                    "p_value_binomtest": 4.5140965972544436e-07,
                    "effect_size_binomtest": 1.3592692844275622,
                },
                "ENSG00000118363": {
                    "p_value_binomtest": 1.3670165462597359e-101,
                    "effect_size_binomtest": 1.2869060309287068,
                },
                "ENSG00000051108": {
                    "p_value_binomtest": 1.7892111295524904e-104,
                    "effect_size_binomtest": 1.1517488666121325,
                },
            }
            self.assertDictEqual(result, expected)
