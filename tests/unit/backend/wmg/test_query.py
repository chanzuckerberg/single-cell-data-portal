import unittest
from typing import NamedTuple

from backend.wmg.api.v1 import get_dot_plot_data, agg_cell_type_counts, agg_tissue_counts
from backend.wmg.data.query import WmgQueryCriteria, WmgQuery
from backend.wmg.data.schemas.cube_schema import cube_non_indexed_dims
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    create_temp_wmg_snapshot,
    all_ones_expression_summary_values,
    all_tens_cell_counts_values,
    all_X_cell_counts_values,
)


# TODO: Test build_* methods separately in test_v1.py.  This package's unit tests need only test the raw results of
#  WmgQuery methods


class QueryTest(unittest.TestCase):
    def test__query_with_no_genes__returns_empty_result(self):
        criteria = WmgQueryCriteria(
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
        )

        dim_size = 1
        with create_temp_wmg_snapshot(
            dim_size=dim_size, expression_summary_vals_fn=all_ones_expression_summary_values
        ) as snapshot:
            query = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(query.expression_summary(criteria), query.cell_counts(criteria))

        expected = {
            "cell_type_ontology_term_id": {},
            "gene_ontology_term_id": {},
            "n_cells": {},
            "n_cells_cell_type": {},
            "n_cells_tissue": {},
            "nnz": {},
            "sum": {},
            "tissue_ontology_term_id": {},
        }
        self.assertEqual(expected, result.to_dict())

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
            query = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(query.expression_summary(criteria), query.cell_counts(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 1)
        assert expected_cell_count_per_cell_type == 729

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
        ]

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
            query = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(query.expression_summary(criteria), query.cell_counts(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 1)
        assert expected_cell_count_per_cell_type == 729

        expected_cell_count_per_tissue = 10 * (dim_size ** len(cube_non_indexed_dims))
        assert expected_cell_count_per_tissue == 21870

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells": 729,
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
        ]

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
            query = WmgQuery(snapshot)
            result = agg_cell_type_counts(query.cell_counts(criteria))

        expected_n_combinations = dim_size ** (len(cube_non_indexed_dims) - 1)
        assert expected_n_combinations == 729

        # after aggregating, we will get three tissues, and three cell types per tissue,
        # with 729 * expected_count total cells
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
            query = WmgQuery(snapshot)
            result = agg_tissue_counts(query.cell_counts(criteria))

        expected_n_combinations = dim_size ** (len(cube_non_indexed_dims) - 1)
        assert expected_n_combinations == 729

        # after aggregating, we will get three tissues,
        # with 729 * expected_count * (# cell types per tissue = 3) total cells
        expected = [{"n_cells_tissue": expected_n_combinations * expected_count * dim_size}] * len(
            criteria.tissue_ontology_term_ids
        )
        self.assertEqual(expected, result.to_dict("records"))

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.maxDiff = None

    def test__query_non_indexed_dim_single_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            dataset_ids=["dataset_id_1"],  # <-- non-indexed dim, single-valued
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            query = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(query.expression_summary(criteria), query.cell_counts(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 2)
        assert expected_cell_count_per_cell_type == 243

        expected_cell_count_per_tissue = 10 * (dim_size ** (len(cube_non_indexed_dims) - 1))
        assert expected_cell_count_per_tissue == 7290

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells": 243,
                "n_cells_cell_type": 2430,
                "n_cells_tissue": 7290,
                "nnz": 243,
                "sum": 243.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells": 243,
                "n_cells_cell_type": 2430,
                "n_cells_tissue": 7290,
                "nnz": 243,
                "sum": 243.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells": 243,
                "n_cells_cell_type": 2430,
                "n_cells_tissue": 7290,
                "nnz": 243,
                "sum": 243.0,
            },
        ]

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
            dataset_ids=["dataset_id_1", "dataset_id_0"],  # <-- non-indexed dim, multi-valued
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            query = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(query.expression_summary(criteria), query.cell_counts(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 2) * 2
        assert expected_cell_count_per_cell_type == 486

        expected_cell_count_per_tissue = 10 * (dim_size ** (len(cube_non_indexed_dims) - 1) * 2)
        assert expected_cell_count_per_tissue == 14580

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells": 486,
                "n_cells_cell_type": 4860,
                "n_cells_tissue": 14580,
                "nnz": 486,
                "sum": 486.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells": 486,
                "n_cells_cell_type": 4860,
                "n_cells_tissue": 14580,
                "nnz": 486,
                "sum": 486.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells": 486,
                "n_cells_cell_type": 4860,
                "n_cells_tissue": 14580,
                "nnz": 486,
                "sum": 486.0,
            },
        ]

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
            ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_1"],  # <-- non-indexed dim, single-valued
            dataset_ids=["dataset_id_1", "dataset_id_0"],  # <-- non-indexed dim, multi-valued
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            query = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(query.expression_summary(criteria), query.cell_counts(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 3) * 1 * 2
        assert expected_cell_count_per_cell_type == 162

        expected_cell_count_per_tissue = 10 * (dim_size ** (len(cube_non_indexed_dims) - 2) * 1 * 2)
        assert expected_cell_count_per_tissue == 4860

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells": 162,
                "n_cells_cell_type": 1620,
                "n_cells_tissue": 4860,
                "nnz": 162,
                "sum": 162.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells": 162,
                "n_cells_cell_type": 1620,
                "n_cells_tissue": 4860,
                "nnz": 162,
                "sum": 162.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells": 162,
                "n_cells_cell_type": 1620,
                "n_cells_tissue": 4860,
                "nnz": 162,
                "sum": 162.0,
            },
        ]

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


class QueryPrimaryFilterDimensionsTest(unittest.TestCase):
    def test__single_dimension__returns_all_dimension_and_terms(self):
        dim_size = 3
        with create_temp_wmg_snapshot(dim_size=dim_size) as snapshot:
            result = WmgQuery(snapshot).list_primary_filter_dimension_term_ids("gene_ontology_term_id")
            self.assertEquals(["gene_ontology_term_id_0", "gene_ontology_term_id_1", "gene_ontology_term_id_2"], result)

    def test__multiple_dimensions__returns_all_dimensions_and_terms_as_tuples(self):
        dim_size = 3

        def exclude_one_gene_per_organism(logical_coord: NamedTuple) -> bool:
            # HACK: method called during building of both "expr summary" and "cell count" cubes, but the latter does not
            # include gene_ontology_term_id
            if "gene_ontology_term_id" not in logical_coord._fields:
                return False
            return logical_coord.gene_ontology_term_id == logical_coord.organism_ontology_term_id.replace(
                "organism", "gene"
            )

        with create_temp_wmg_snapshot(
            dim_size=dim_size, exclude_logical_coord_fn=exclude_one_gene_per_organism
        ) as snapshot:
            result = WmgQuery(snapshot).list_grouped_primary_filter_dimensions_term_ids(
                "gene_ontology_term_id", "organism_ontology_term_id"
            )
            self.assertEquals(
                {
                    "organism_ontology_term_id_0": ["gene_ontology_term_id_1", "gene_ontology_term_id_2"],
                    "organism_ontology_term_id_1": ["gene_ontology_term_id_0", "gene_ontology_term_id_2"],
                    "organism_ontology_term_id_2": ["gene_ontology_term_id_0", "gene_ontology_term_id_1"],
                },
                result,
            )
