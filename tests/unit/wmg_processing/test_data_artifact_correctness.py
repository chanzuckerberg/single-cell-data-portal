import unittest
from typing import Dict

from backend.wmg.data.query import WmgQueryCriteria, WmgSnapshot
from backend.wmg.data.utils import find_all_dim_option_values, find_dim_option_values
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot

TEST_SNAPSHOT = "realistic-test-snapshot"
# disease ontolgoy term id for normal cells
NORMAL_CELL_DISEASE_ONTOLOGY_TERM_ID = "PATO:0000461"


def find_dim_option_values_pandas(criteria: Dict, snapshot: WmgSnapshot, dimension: str) -> set:
    """Find values for the specified dimension that satisfy the given filtering criteria,
    ignoring any criteria specified for the given dimension."""
    filter_options_criteria = criteria.copy(update={dimension + "s": []}, deep=True)
    cell_counts_df = snapshot.cell_counts_cube.df[:]

    for key, vals in filter_options_criteria:
        key = key[:-1] if key[-1] == "s" else key
        vals = [vals] if isinstance(vals, str) else vals
        if key in cell_counts_df.columns and len(vals) > 0:
            cell_counts_df = cell_counts_df[cell_counts_df[key].isin(vals)]

    filter_dims = list(cell_counts_df.groupby(dimension).groups.keys())
    return filter_dims


class DataArtifactCorrectness(unittest.TestCase):
    # creates filter dim options by querying the cell counts cube and checks to make sure that the result
    # is equal to the current approach of creating filter dim options via the precomputed filter relationships
    def test__filter_dim_options_correctness_all_filter_values(self):
        # this tests filter options equality for the entire data corpus
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            cell_counts = snapshot.cell_counts_cube.df[:]
            cell_counts = cell_counts.select_dtypes(exclude="number")
            cell_counts = cell_counts[cell_counts["organism_ontology_term_id"] == "NCBITaxon:9606"]
            expected_filter_options = {
                dimension: list(set(cell_counts[dimension]))
                for dimension in cell_counts
                if dimension != "organism_ontology_term_id"
            }
            test_filter_options = {
                dimension: find_all_dim_option_values(snapshot, "NCBITaxon:9606", dimension)
                for dimension in cell_counts
                if dimension != "organism_ontology_term_id"
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
                dims_expected[dim] = find_dim_option_values_pandas(criteria, snapshot, dim)
                dims[dim].sort()
                dims_expected[dim].sort()

            self.assertDictEqual(dims, dims_expected)
