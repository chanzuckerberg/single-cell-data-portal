import json
import unittest
from math import log
from unittest.mock import patch

from backend.de.server.app import app
from tests.unit.backend.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class DeAPIV1Tests(unittest.TestCase):
    def setUp(self):
        super().setUp()
        with EnvironmentSetup(dict(APP_NAME="corpora-api-de")):
            self.app = app.test_client(use_cookies=False)

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.maxDiff = None

    @patch("backend.de.api.v1.load_snapshot")
    def test__differentialExpression_returns_expected_results(self, load_snapshot):
        expected_effect_size_sums = [
            [-334, 648, 71, 2654, 696, 2763, 2021],
            [-334, 648, 71, 2654, 696, 2763, 0],
            [-334, 648, 71, 2654, 696, 2763, 2836],
        ]
        expected_log_fold_change_sums = [
            [-54, 102, 66, 491, 199, 532, 483],
            [-54, 102, 66, 491, 199, 532, 0],
            [-54, 102, 66, 491, 199, 532, 500],
        ]
        expected_log_p_value_sums = [
            [348318, 133077, 195787, 125959, 134096, 69440, 123317],
            [348318, 133077, 195787, 125959, 134096, 69440, 0],
            [348318, 133077, 195787, 125959, 134096, 69440, 123599],
        ]
        expected_n_overlap = [
            [0, 0, 0, 0, 0, 0, 37],
            [0, 0, 0, 0, 0, 0, 37],
            [0, 0, 0, 0, 0, 0, 37],
        ]

        for test_index, excludeOverlappingCells in enumerate(["retainBoth", "excludeOne", "excludeTwo"]):
            test_cases = [
                {
                    "queryGroup1Filters": {
                        "organism_ontology_term_id": "NCBITaxon:9606",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": ["CL:0000066"],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": [],
                    },
                    "queryGroup2Filters": {
                        "organism_ontology_term_id": "NCBITaxon:9606",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": ["CL:0000084"],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": [],
                    },
                },
                {
                    "queryGroup1Filters": {
                        "organism_ontology_term_id": "NCBITaxon:9606",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": [],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0014"],
                        "sex_ontology_term_ids": [],
                    },
                    "queryGroup2Filters": {
                        "organism_ontology_term_id": "NCBITaxon:9606",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": [],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0005"],
                        "sex_ontology_term_ids": [],
                    },
                },
                {
                    "queryGroup1Filters": {
                        "organism_ontology_term_id": "NCBITaxon:9606",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": ["CL:0000066"],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0014"],
                        "sex_ontology_term_ids": [],
                    },
                    "queryGroup2Filters": {
                        "organism_ontology_term_id": "NCBITaxon:9606",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": ["CL:0000084"],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0005"],
                        "sex_ontology_term_ids": [],
                    },
                },
                {
                    "queryGroup1Filters": {
                        "organism_ontology_term_id": "NCBITaxon:10090",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": ["CL:0000115"],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": [],
                    },
                    "queryGroup2Filters": {
                        "organism_ontology_term_id": "NCBITaxon:10090",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": ["CL:0000169"],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": [],
                    },
                },
                {
                    "queryGroup1Filters": {
                        "organism_ontology_term_id": "NCBITaxon:10090",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": [],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": ["PATO:0000383"],
                    },
                    "queryGroup2Filters": {
                        "organism_ontology_term_id": "NCBITaxon:10090",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": [],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": ["PATO:0000384"],
                    },
                },
                {
                    "queryGroup1Filters": {
                        "organism_ontology_term_id": "NCBITaxon:10090",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": ["CL:0000115"],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": ["PATO:0000383"],
                    },
                    "queryGroup2Filters": {
                        "organism_ontology_term_id": "NCBITaxon:10090",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": ["CL:0000169"],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": ["PATO:0000384"],
                    },
                },
                {
                    "queryGroup1Filters": {
                        "organism_ontology_term_id": "NCBITaxon:10090",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": ["CL:0000115"],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": ["PATO:0000383"],
                    },
                    "queryGroup2Filters": {
                        "organism_ontology_term_id": "NCBITaxon:10090",
                        "tissue_ontology_term_ids": [],
                        "cell_type_ontology_term_ids": [],
                        "publication_citations": [],
                        "disease_ontology_term_ids": [],
                        "self_reported_ethnicity_ontology_term_ids": [],
                        "sex_ontology_term_ids": ["PATO:0000383"],
                    },
                },
            ]

            with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
                for i, test_case in enumerate(test_cases):
                    test_case["excludeOverlappingCells"] = excludeOverlappingCells

                    with self.subTest(test_case=test_case):
                        load_snapshot.return_value = snapshot
                        response = self.app.post(
                            "/de/v1/differentialExpression",
                            headers={"Content-Type": "application/json"},
                            data=json.dumps(test_case),
                        )
                        self.assertEqual(response.status_code, 200)
                        result = json.loads(response.data)
                        effect_size_sum = round(
                            sum([i["effect_size"] for i in result["differentialExpressionResults"]])
                        )
                        log_p_value_sum = round(
                            sum([-log(i["adjusted_p_value"] + 1e-300) for i in result["differentialExpressionResults"]])
                        )
                        log_fold_change_sum = round(
                            sum([i["log_fold_change"] for i in result["differentialExpressionResults"]])
                        )

                        self.assertEqual(effect_size_sum, expected_effect_size_sums[test_index][i])
                        self.assertEqual(log_p_value_sum, expected_log_p_value_sums[test_index][i])
                        self.assertEqual(log_fold_change_sum, expected_log_fold_change_sums[test_index][i])
                        self.assertEqual(result["n_overlap"], expected_n_overlap[test_index][i])
