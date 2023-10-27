"""This module tests the helper functions used in `backend.wmg.api.v2.py`.
"""

import unittest

from parameterized import parameterized

from backend.wmg.api.v2 import build_ontology_term_id_label_mapping


class ExpressionDotPlotTest(unittest.TestCase):
    def _ontology_term_id_to_label_map_test_cases():
        test_cases = [
            (
                "disease_ontology_term_id",
                ["MONDO:0100096", "MONDO:0000001"],
                [{"MONDO:0100096": "COVID-19"}, {"MONDO:0000001": "disease"}],
            ),
            (
                "tissue_ontology_term_id",
                ["UBERON:0002539", "UBERON:0000995 (organoid)"],
                [{"UBERON:0002539": "pharyngeal arch"}, {"UBERON:0000995 (organoid)": "uterus (organoid)"}],
            ),
            (
                "sex_ontology_term_id",
                ["PATO:0000384", "PATO:0000383", "unknown"],
                [{"PATO:0000384": "male"}, {"PATO:0000383": "female"}, {"unknown": None}],
            ),
            (
                "self_reported_ethnicity_ontology_term_id",
                [
                    # Schema-4 ontology term ids
                    "HANCESTRO:0005,HANCESTRO:0014",
                    # Schema-3 ontology term ids
                    "HANCESTRO:0003",
                    "unknown",
                    "multiethnic",
                    "na",
                ],
                [
                    # Schema-4 values
                    {"HANCESTRO:0005,HANCESTRO:0014": "European,Hispanic or Latin American"},
                    # Schema-3 values
                    {"HANCESTRO:0003": "country"},
                    {"unknown": None},
                    {"multiethnic": None},
                    {"na": None},
                ],
            ),
            ("development_stage_ontology_term_id", ["HsapDv:0000000"], [{"HsapDv:0000000": "human life cycle stage"}]),
            (
                "cell_type_ontology_term_id",
                ["CL:0000000", "CL:0000082 (cell culture)"],
                [{"CL:0000000": "cell"}, {"CL:0000082 (cell culture)": "epithelial cell of lung (cell culture)"}],
            ),
        ]

        return test_cases

    @parameterized.expand(_ontology_term_id_to_label_map_test_cases)
    def test__build_ontology_term_id_label_mapping(self, dimension, input_ontology_term_id_list, expected_output_list):
        actual_term_id_to_label = build_ontology_term_id_label_mapping(
            ontology_term_ids=input_ontology_term_id_list, dim_name=dimension
        )

        self.assertEqual(actual_term_id_to_label, expected_output_list)
