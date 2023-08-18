"""This module tests the functions that build the gene expression and cell counts data structures.

In detail, this module tests the public and private functions defined in 
`backend.wmg.api.common.expression_dotplot_data` module.
"""

import unittest

from backend.wmg.api.common.expression_dotplot import agg_cell_type_counts, agg_tissue_counts, get_dot_plot_data
from backend.wmg.api.wmg_api_config import (
    READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
    READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
)
from backend.wmg.data.query import WmgCubeQueryParams, WmgQuery, WmgQueryCriteria
from tests.unit.backend.wmg.fixtures.test_cube_schema import expression_summary_non_indexed_dims
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    all_ones_expression_summary_values,
    all_tens_cell_counts_values,
    all_X_cell_counts_values,
    create_temp_wmg_snapshot,
)

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


def generate_expected_dot_plot_data_with_pandas(snapshot, criteria):
    """
    Build expected query results from a pandas dataframe.
    """
    expression_summary = snapshot.expression_summary_cube.df[:]
    cell_counts = snapshot.cell_counts_cube.df[:]

    criteria_es = criteria.dict()
    criteria_cc = criteria.copy(exclude={"gene_ontology_term_ids"}).dict()
    expression_summary = _filter_dataframe(expression_summary, criteria_es)
    cell_counts = _filter_dataframe(cell_counts, criteria_cc)
    cell_counts.rename(columns={"n_cells": "n_total_cells"}, inplace=True)
    expected, _ = get_dot_plot_data(expression_summary, cell_counts)
    return sorted(
        expected.to_dict("records"),
        key=lambda r: (
            r["gene_ontology_term_id"],
            r["tissue_ontology_term_id"],
            r["cell_type_ontology_term_id"],
        ),
    )


class ExpressionDotPlotTest(unittest.TestCase):
    def setUp(self):
        self.cube_query_params = WmgCubeQueryParams(
            cube_query_valid_attrs=READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
            cube_query_valid_dims=READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
        )

    def test__query_all_indexed_dims_single_value__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_1",
            tissue_ontology_term_ids=["tissue_ontology_term_id_2"],
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

            expected = generate_expected_dot_plot_data_with_pandas(snapshot, criteria)

            self.assertEqual(
                expected,
                sorted(
                    result.to_dict("records"),
                    key=lambda r: (
                        r["gene_ontology_term_id"],
                        r["tissue_ontology_term_id"],
                        r["cell_type_ontology_term_id"],
                    ),
                ),
            )

    def test__query_all_indexed_dims_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0", "gene_ontology_term_id_2"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"],
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

            expected = generate_expected_dot_plot_data_with_pandas(snapshot, criteria)

            self.assertEqual(
                expected,
                sorted(
                    result.to_dict("records"),
                    key=lambda r: (
                        r["gene_ontology_term_id"],
                        r["tissue_ontology_term_id"],
                        r["cell_type_ontology_term_id"],
                    ),
                ),
            )

    def test__query_non_indexed_dim_single_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            development_stage_ontology_term_ids=[
                "development_stage_ontology_term_id_1"
            ],  # <-- non-indexed dim, single-valued
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

            expected = generate_expected_dot_plot_data_with_pandas(snapshot, criteria)

            self.assertEqual(
                expected,
                sorted(
                    result.to_dict("records"),
                    key=lambda r: (
                        r["gene_ontology_term_id"],
                        r["tissue_ontology_term_id"],
                        r["cell_type_ontology_term_id"],
                    ),
                ),
            )

    def test__query_non_indexed_dim_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            development_stage_ontology_term_ids=[
                "development_stage_ontology_term_id_1",
                "development_stage_ontology_term_id_0",
            ],  # <-- non-indexed dim, multi-valued
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

            expected = generate_expected_dot_plot_data_with_pandas(snapshot, criteria)

            self.assertEqual(
                expected,
                sorted(
                    result.to_dict("records"),
                    key=lambda r: (
                        r["gene_ontology_term_id"],
                        r["tissue_ontology_term_id"],
                        r["cell_type_ontology_term_id"],
                    ),
                ),
            )

    def test__query_non_indexed_dim_single_and_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            self_reported_ethnicity_ontology_term_ids=[
                "self_reported_ethnicity_ontology_term_id_1"
            ],  # <-- non-indexed dim, single-valued
            development_stage_ontology_term_ids=[
                "development_stage_ontology_term_id_1",
                "development_stage_ontology_term_id_0",
            ],  # <-- non-indexed dim, multi-valued
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

            expected = generate_expected_dot_plot_data_with_pandas(snapshot, criteria)
            self.assertEqual(
                expected,
                sorted(
                    result.to_dict("records"),
                    key=lambda r: (
                        r["gene_ontology_term_id"],
                        r["tissue_ontology_term_id"],
                        r["cell_type_ontology_term_id"],
                    ),
                ),
            )

    def test__query_agg_cell_type_counts__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_1",
            tissue_ontology_term_ids=[
                "tissue_ontology_term_id_0",
                "tissue_ontology_term_id_1",
                "tissue_ontology_term_id_2",
            ],
        )

        dim_size = 3
        expected_count = 42
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=lambda coords: all_X_cell_counts_values(coords, expected_count),
        ) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
            result = agg_cell_type_counts(q.cell_counts(criteria))

        not_used_cube_indexed_dims = [0 if criteria.dict()[dim_name] else 1 for dim_name in ALL_INDEXED_DIMS_FOR_QUERY]
        expected_n_combinations = dim_size ** (
            len(expression_summary_non_indexed_dims) + sum(not_used_cube_indexed_dims) - 1
        )

        expected = (
            [{"n_cells_cell_type": expected_n_combinations * expected_count}]
            * len(criteria.tissue_ontology_term_ids)
            * dim_size
        )
        self.assertEqual(expected, result.to_dict("records"))

    def test__query_agg_tissue_counts__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_1",
            tissue_ontology_term_ids=[
                "tissue_ontology_term_id_0",
                "tissue_ontology_term_id_1",
                "tissue_ontology_term_id_2",
            ],
        )

        dim_size = 3
        expected_count = 42
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=lambda coords: all_X_cell_counts_values(coords, expected_count),
        ) as snapshot:
            q = WmgQuery(snapshot, self.cube_query_params)
            result = agg_tissue_counts(q.cell_counts(criteria))

        not_used_cube_indexed_dims = [0 if criteria.dict()[dim_name] else 1 for dim_name in ALL_INDEXED_DIMS_FOR_QUERY]
        expected_n_combinations = dim_size ** (
            len(expression_summary_non_indexed_dims) + sum(not_used_cube_indexed_dims) - 1
        )

        expected = [
            {
                "n_cells_tissue": expected_n_combinations * expected_count * dim_size,
            }
        ] * len(criteria.tissue_ontology_term_ids)
        self.assertEqual(expected, result.to_dict("records"))
