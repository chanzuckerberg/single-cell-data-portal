import unittest
from typing import NamedTuple

from backend.wmg.data.query import (
    FmgQueryCriteria,
    MarkerGeneQueryCriteria,
    WmgQuery,
    WmgQueryCriteria,
    retrieve_top_n_markers,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    create_temp_wmg_snapshot,
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"

ALL_INDEXED_DIMS_FOR_QUERY = [
    "gene_ontology_term_ids",
    "tissue_ontology_term_ids",
    "organism_ontology_term_id",
]

# TODO: Test build_* methods separately in test_v1.py.  This package's unit tests need only test the raw results of
#  query methods


def _filter_dataframe(dataframe, criteria):
    for key in criteria:
        attrs = [criteria[key]] if not isinstance(criteria[key], list) else criteria[key]
        if len(attrs) > 0:
            depluralized_key = key[:-1] if key[-1] == "s" else key
            dataframe = dataframe[dataframe[depluralized_key].isin(attrs)]
    return dataframe


def generate_expected_marker_gene_data_with_pandas(snapshot, criteria, statistical_test, num_markers):
    """
    Build expected query results from a pandas dataframe.
    """
    marker_genes = snapshot.marker_genes_cube.df[:]

    criteria_mg = criteria.dict()
    marker_genes = _filter_dataframe(marker_genes, criteria_mg)
    expected = retrieve_top_n_markers(marker_genes, statistical_test, num_markers)
    return expected


class QueryTest(unittest.TestCase):
    def test__query_marker_genes_cube__returns_correct_top_10_markers(self):
        criteria = MarkerGeneQueryCriteria(
            tissue_ontology_term_id="UBERON:0002048",
            cell_type_ontology_term_id="CL:0000786",
            organism_ontology_term_id="NCBITaxon:9606",
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            result = q.marker_genes(criteria)
            marker_genes = retrieve_top_n_markers(result, "ttest", 10)
            expected = generate_expected_marker_gene_data_with_pandas(snapshot, criteria, "ttest", 10)
            self.assertEqual(marker_genes, expected)

    def test__query_marker_genes_cube__returns_correct_all_markers(self):
        criteria = MarkerGeneQueryCriteria(
            tissue_ontology_term_id="UBERON:0002048",
            cell_type_ontology_term_id="CL:0000786",
            organism_ontology_term_id="NCBITaxon:9606",
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            result = q.marker_genes(criteria)
            marker_genes = retrieve_top_n_markers(result, "ttest", 0)
            expected = generate_expected_marker_gene_data_with_pandas(snapshot, criteria, "ttest", 0)
            self.assertEqual(marker_genes, expected)

    def test__query_expression_summary_fmg_cube__returns_correct_results(self):
        criteria = FmgQueryCriteria(
            gene_ontology_term_ids=["ENSG00000238042", "ENSG00000168028"],
            organism_ontology_term_id="NCBITaxon:9606",
            tissue_ontology_term_ids=["UBERON:0002048"],
            cell_type_ontology_term_ids=["CL:0000786"],
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            query_result = q.expression_summary_fmg(criteria)
            query_sum = list(query_result[["sum", "sqsum", "nnz", "nnz_thr"]].sum())
            expected = [120129.5703125, 343018.5625, 49966.0, 48978.0]
            [self.assertEqual(round(query_sum[i]), round(expected[i])) for i in range(len(query_sum))]

    def test__query_expression_summary_default_cube__returns_correct_results(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["ENSG00000238042", "ENSG00000168028"],
            organism_ontology_term_id="NCBITaxon:9606",
            tissue_ontology_term_ids=["UBERON:0002048"],
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            query_result = q.expression_summary_default(criteria)
            query_sum = list(query_result[["sum", "nnz", "sqsum"]].sum())
            expected = [804599.0, 370972.0, 1822108.0]
            [self.assertEqual(round(query_sum[i]), round(expected[i])) for i in range(len(query_sum))]


class QueryPrimaryFilterDimensionsTest(unittest.TestCase):
    def test__single_dimension__returns_all_dimension_and_terms(self):
        dim_size = 3
        with create_temp_wmg_snapshot(dim_size=dim_size) as snapshot:
            q = WmgQuery(snapshot)
            result = q.list_primary_filter_dimension_term_ids("tissue_ontology_term_id")
            self.assertEquals(
                ["tissue_ontology_term_id_0", "tissue_ontology_term_id_1", "tissue_ontology_term_id_2"], result
            )

    def test__multiple_dimensions__returns_all_dimensions_and_terms_as_tuples(self):
        dim_size = 3

        def exclude_one_tissue_per_organism(logical_coord: NamedTuple) -> bool:
            return logical_coord.tissue_ontology_term_id == logical_coord.organism_ontology_term_id.replace(
                "organism", "tissue"
            )

        with create_temp_wmg_snapshot(
            dim_size=dim_size, exclude_logical_coord_fn=exclude_one_tissue_per_organism
        ) as snapshot:
            q = WmgQuery(snapshot)
            result = q.list_grouped_primary_filter_dimensions_term_ids(
                "tissue_ontology_term_id", "organism_ontology_term_id"
            )
            self.assertEquals(
                {
                    "organism_ontology_term_id_0": ["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"],
                    "organism_ontology_term_id_1": ["tissue_ontology_term_id_0", "tissue_ontology_term_id_2"],
                    "organism_ontology_term_id_2": ["tissue_ontology_term_id_0", "tissue_ontology_term_id_1"],
                },
                result,
            )
