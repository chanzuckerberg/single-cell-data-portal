import unittest
from typing import Dict

from pandas.testing import assert_series_equal

from backend.wmg.data.query import WmgQuery, WmgQueryCriteria, WmgSnapshot
from backend.wmg.data.schemas.cube_schema_default import (
    expression_summary_logical_dims as expression_summary_default_logical_dims,
)
from backend.wmg.data.utils import find_all_dim_option_values, find_dim_option_values
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot

TEST_SNAPSHOT = "realistic-test-snapshot"
# disease ontolgoy term id for normal cells
NORMAL_CELL_DISEASE_ONTOLOGY_TERM_ID = "PATO:0000461"


def find_dim_option_values_tiledb(criteria: Dict, snapshot: WmgSnapshot, dimension: str) -> set:
    """Find values for the specified dimension that satisfy the given filtering criteria,
    ignoring any criteria specified for the given dimension."""
    filter_options_criteria = criteria.copy(update={dimension + "s": []}, deep=True)
    # todo can we query cell_counts for a performance gain?
    q = WmgQuery(snapshot)
    query_result = q.cell_counts(filter_options_criteria)
    filter_dims = list(query_result.groupby(dimension).groups.keys())
    return filter_dims


class DataArtifactCorrectness(unittest.TestCase):
    # creates the expression summary default cube by group-by and summing the expression summary cube
    # and checks to makes sure that the result is the same as the existing fixture (expression_summary_default)
    def test__expression_summary_default_cube_correctness(self):
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            expression_summary_df = snapshot.expression_summary_cube.df[:]
            expression_summary_default_df = snapshot.expression_summary_default_cube.df[:]

            # generate expression summary default cube by group-by and summing the expression summary cube
            test_expression_summary_default_df = (
                expression_summary_df.groupby(expression_summary_default_logical_dims)
                .sum(numeric_only=True)
                .reset_index()
            )

            # check that the two cubes are equal
            df1 = expression_summary_default_df.sort_values(expression_summary_default_logical_dims).reset_index(
                drop=True
            )
            df2 = test_expression_summary_default_df.sort_values(expression_summary_default_logical_dims).reset_index(
                drop=True
            )
            [assert_series_equal(df1[col], df2[col], rtol=1e-2) for col in df1]

    # creates filter dim options by querying the cell counts cube and checks to make sure that the result
    # is equal to the current approach of creating filter dim options via the precomputed filter relationships
    def test__filter_dim_options_correctness_all_filter_values(self):
        # this tests filter options equality for the entire data corpus
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            cell_counts = snapshot.cell_counts_cube.df[:]
            cell_counts = cell_counts.select_dtypes(exclude="number")
            expected_filter_options = {dimension: list(set(cell_counts[dimension])) for dimension in cell_counts}
            test_filter_options = {
                dimension: find_all_dim_option_values(snapshot, dimension) for dimension in cell_counts
            }

            for dimension in test_filter_options:
                test_filter_options[dimension].sort()
                expected_filter_options[dimension].sort()

            # check that the two dictionaries are equal
            self.assertDictEqual(test_filter_options, expected_filter_options)

    def test__filter_dim_options_correctness_with_criteria(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["ENSG00000182149"],
            organism_ontology_term_id="NCBITaxon:9606",
            disease_ontology_term_ids=[NORMAL_CELL_DISEASE_ONTOLOGY_TERM_ID],
            tissue_ontology_term_ids=["UBERON:0002048"],
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            dims = {
                "dataset_id": "",
                "disease_ontology_term_id": "",
                "sex_ontology_term_id": "",
                "development_stage_ontology_term_id": "",
                "self_reported_ethnicity_ontology_term_id": "",
                "tissue_ontology_term_id": "",
            }
            dims_expected = dims.copy()
            for dim in dims:
                dims[dim] = find_dim_option_values(criteria, snapshot, dim)
                dims_expected[dim] = find_dim_option_values_tiledb(criteria, snapshot, dim)
                dims[dim].sort()
                dims_expected[dim].sort()

            self.assertDictEqual(dims, dims_expected)
