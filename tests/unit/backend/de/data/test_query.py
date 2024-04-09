import unittest

from backend.wmg.data.query import DeQueryCriteria, _select_cube_with_best_discriminatory_power
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot

TEST_SNAPSHOT = "realistic-test-snapshot"


class DeQueryTests(unittest.TestCase):
    def test__correct_cube_selected_based_on_query(self):
        test_cases = [
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": ["CL:0000066"],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": [],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": ["CL:0000084"],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": [],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": [],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0014"],
                "sex_ontology_term_ids": [],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": [],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0005"],
                "sex_ontology_term_ids": [],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": ["CL:0000066"],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0014"],
                "sex_ontology_term_ids": [],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": ["CL:0000084"],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0005"],
                "sex_ontology_term_ids": [],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:10090",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": ["CL:0000115"],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": [],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:10090",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": ["CL:0000169"],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": [],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:10090",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": [],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": ["PATO:0000383"],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:10090",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": [],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": ["PATO:0000384"],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:10090",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": ["CL:0000115"],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": ["PATO:0000383"],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:10090",
                "tissue_ontology_term_ids": [],
                "cell_type_ontology_term_ids": ["CL:0000169"],
                "publication_citations": [],
                "disease_ontology_term_ids": [],
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": ["PATO:0000384"],
            },
        ]
        expected_cube_keys = [
            "default",
            "default",
            "self_reported_ethnicity_ontology_term_id",
            "self_reported_ethnicity_ontology_term_id",
            "self_reported_ethnicity_ontology_term_id",
            "self_reported_ethnicity_ontology_term_id",
            "default",
            "default",
            "sex_ontology_term_id",
            "sex_ontology_term_id",
            "sex_ontology_term_id",
            "sex_ontology_term_id",
        ]
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            for i, test_case in enumerate(test_cases):
                with self.subTest(test_case=test_case):
                    criteria = DeQueryCriteria(**test_case)
                    cube = _select_cube_with_best_discriminatory_power(snapshot, criteria)
                    cube_key = cube.uri.split("__")[-1]
                    self.assertEqual(cube_key, expected_cube_keys[i])
