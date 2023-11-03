"""This module tests the helper functions used in `backend.wmg.api.v2.py`.
"""

import pytest

from backend.wmg.api.v2 import sanitize_api_query_dict


@pytest.mark.parametrize(
    "input_query,expected_sanitized_query",
    [
        (
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0008"],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0008"],
            },
        ),
        (
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0008", "HANCESTRO:0008,HANCESTOR:0021"],
            },
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0008"],
            },
        ),
        (
            {
                "organism_ontology_term_id": "NCBITaxon:9606",
                "self_reported_ethnicity_ontology_term_ids": ["HANCESTRO:0008,HANCESTOR:0021"],
            },
            {"organism_ontology_term_id": "NCBITaxon:9606", "self_reported_ethnicity_ontology_term_ids": []},
        ),
        (
            {"organism_ontology_term_id": "NCBITaxon:9606", "self_reported_ethnicity_ontology_term_ids": []},
            {"organism_ontology_term_id": "NCBITaxon:9606", "self_reported_ethnicity_ontology_term_ids": []},
        ),
        ({"organism_ontology_term_id": "NCBITaxon:9606"}, {"organism_ontology_term_id": "NCBITaxon:9606"}),
    ],
)
def test_sanitize_api_query_dict(input_query, expected_sanitized_query):
    # NOTE: `sanitize_api_query_dict()` mutates the function argument
    sanitize_api_query_dict(input_query)

    assert input_query == expected_sanitized_query
