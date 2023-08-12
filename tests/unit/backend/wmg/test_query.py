import unittest
from typing import NamedTuple

from backend.wmg.api.wmg_api_config import (
    READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
    READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
)
from backend.wmg.data.query import (
    FmgQueryCriteria,
    MarkerGeneQueryCriteria,
    WmgCubeQueryParams,
    WmgQuery,
    WmgQueryCriteria,
    retrieve_top_n_markers,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import create_temp_wmg_snapshot, load_realistic_test_snapshot

TEST_SNAPSHOT = "realistic-test-snapshot"

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


class WmgCubeQueryParamsTest(unittest.TestCase):
    def test__cube_query_params_for_expression_summary_cube(self):
        cube_query_valid_attrs = ["gene_ontology_term_id", "cell_type_ontology_term_id"]
        cube_query_valid_dims = ["gene_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"]

        cube_query_params = WmgCubeQueryParams(
            cube_query_valid_attrs=cube_query_valid_attrs,
            cube_query_valid_dims=cube_query_valid_dims,
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            expected_cube_indexed_dims = [
                "gene_ontology_term_id",
                "tissue_ontology_term_id",
                "organism_ontology_term_id",
            ]
            actual_cube_indexed_dims = [i.name for i in snapshot.expression_summary_cube.schema.domain]
            self.assertEqual(actual_cube_indexed_dims, expected_cube_indexed_dims)

            expected_cube_nonindexed_attrs = [
                "cell_type_ontology_term_id",
                "tissue_original_ontology_term_id",
                "dataset_id",
                "assay_ontology_term_id",
                "development_stage_ontology_term_id",
                "disease_ontology_term_id",
                "self_reported_ethnicity_ontology_term_id",
                "sex_ontology_term_id",
                "publication_citation",
                "nnz",
                "sum",
            ]
            actual_cube_nonindexed_attrs = [i.name for i in snapshot.expression_summary_cube.schema]
            self.assertEqual(actual_cube_nonindexed_attrs, expected_cube_nonindexed_attrs)

            expected_dim_key_names_to_lookup_query_criteria = [
                "gene_ontology_term_ids",
                "tissue_ontology_term_ids",
                "organism_ontology_term_id",
            ]
            actual_dim_key_names_to_lookup_query_criteria = cube_query_params.get_indexed_dims_to_lookup_query_criteria(
                snapshot.expression_summary_cube, pluralize=True
            )
            self.assertEqual(
                actual_dim_key_names_to_lookup_query_criteria, expected_dim_key_names_to_lookup_query_criteria
            )

            expected_attrs_for_cube_query = ["cell_type_ontology_term_id"]
            actual_attrs_for_cube_query = cube_query_params.get_attrs_for_cube_query(snapshot.expression_summary_cube)
            self.assertEqual(actual_attrs_for_cube_query, expected_attrs_for_cube_query)

            expected_dims_for_cube_query = ["gene_ontology_term_id", "tissue_ontology_term_id"]
            actual_dims_for_cube_query = cube_query_params.get_dims_for_cube_query(snapshot.expression_summary_cube)
            self.assertEqual(actual_dims_for_cube_query, expected_dims_for_cube_query)


class QueryTest(unittest.TestCase):
    def setUp(self):
        self.cube_query_params = WmgCubeQueryParams(
            cube_query_valid_attrs=READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
            cube_query_valid_dims=READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
        )

    def test__query_marker_genes_cube__returns_correct_top_10_markers(self):
        criteria = MarkerGeneQueryCriteria(
            tissue_ontology_term_id="UBERON:0002048",
            cell_type_ontology_term_id="CL:0000786",
            organism_ontology_term_id="NCBITaxon:9606",
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
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
            q = WmgQuery(snapshot, self.cube_query_params)
            result = q.marker_genes(criteria)
            marker_genes = retrieve_top_n_markers(result, "ttest", 0)
            expected = generate_expected_marker_gene_data_with_pandas(snapshot, criteria, "ttest", 0)
            self.assertEqual(marker_genes, expected)

    def test__query_expression_summary_fmg_cube__returns_correct_results(self):
        criteria = FmgQueryCriteria(
            organism_ontology_term_id="NCBITaxon:9606",
            tissue_ontology_term_ids=["UBERON:0002048"],
            cell_type_ontology_term_ids=["CL:0000786"],
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
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
            q = WmgQuery(snapshot, self.cube_query_params)
            query_result = q.expression_summary_default(criteria)
            query_sum = list(query_result[["sum", "nnz", "sqsum"]].sum())
            expected = [804599.0, 370972.0, 1822108.0]
            [self.assertEqual(round(query_sum[i]), round(expected[i])) for i in range(len(query_sum))]


class QueryPrimaryFilterDimensionsTest(unittest.TestCase):
    def setUp(self):
        self.cube_query_params = WmgCubeQueryParams(
            cube_query_valid_attrs=READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
            cube_query_valid_dims=READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
        )

    def test__single_dimension__returns_all_dimension_and_terms(self):
        dim_size = 3
        with create_temp_wmg_snapshot(dim_size=dim_size) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
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
            q = WmgQuery(snapshot, self.cube_query_params)
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
