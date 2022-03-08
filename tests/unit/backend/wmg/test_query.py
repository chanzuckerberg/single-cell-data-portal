import unittest
from typing import Tuple

from backend.wmg.data.query import WmgQueryCriteria, WmgQuery, build_dot_plot_matrix
from backend.wmg.data.schema import cube_non_indexed_dims
from tests.unit.backend.wmg.fixtures.cube import create_temp_cube, all_ones_attr_values


class QueryTest(unittest.TestCase):
    def test__query_with_no_genes__returns_empty_result(self):
        criteria = WmgQueryCriteria(
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = build_dot_plot_matrix(WmgQuery(cube).expression_summary(criteria))

        expected = {
            "n_cells": {},
            "nnz": {},
            "sum": {},
        }

        self.assertEqual(expected, result.to_dict())

    def test__query_all_indexed_dims_single_value__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_1",
            tissue_ontology_term_ids=["tissue_ontology_term_id_2"],
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = build_dot_plot_matrix(WmgQuery(cube).expression_summary(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 1)
        assert expected_cell_count_per_cell_type == 729

        expected = {
            "n_cells": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
            },
            "nnz": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type
                ),
            },
        }
        self.maxDiff = None
        self.assertEqual(expected, result.to_dict())

    def test__query_all_indexed_dims_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0", "gene_ontology_term_id_2"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"],
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = build_dot_plot_matrix(WmgQuery(cube).expression_summary(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 1)
        assert expected_cell_count_per_cell_type == 729

        expected = {
            "n_cells": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
            },
            "nnz": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_1",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_2",
                    "tissue_ontology_term_id_2",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type
                ),
            },
        }

        self.assertEqual(expected, result.to_dict())

    def test__query_non_indexed_dim_single_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            dataset_ids=["dataset_id_1"],  # <-- non-indexed dim, single-valued
        )

        dim_size = 2
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = build_dot_plot_matrix(WmgQuery(cube).expression_summary(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 2)
        assert expected_cell_count_per_cell_type == 32

        expected = {
            "n_cells": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
            },
            "nnz": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type
                ),
            },
        }

        self.assertEqual(expected, result.to_dict())

    def test__query_non_indexed_dim_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            dataset_ids=["dataset_id_1", "dataset_id_0"],  # <-- non-indexed dim, multi-valued
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = build_dot_plot_matrix(WmgQuery(cube).expression_summary(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 2) * 2
        assert expected_cell_count_per_cell_type == 486

        expected = {
            "n_cells": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
            },
            "nnz": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type
                ),
            },
        }

        self.assertEqual(expected, result.to_dict())

    def test__query_non_indexed_dim_single_and_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_1"],  # <-- non-indexed dim, single-valued
            dataset_ids=["dataset_id_1", "dataset_id_0"],  # <-- non-indexed dim, multi-valued
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = build_dot_plot_matrix(WmgQuery(cube).expression_summary(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        expected_cell_count_per_cell_type = dim_size ** (len(cube_non_indexed_dims) - 3) * 1 * 2
        assert expected_cell_count_per_cell_type == 162

        expected = {
            "n_cells": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
            },
            "nnz": {
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_0",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_1",
                ): expected_cell_count_per_cell_type,
                (
                    "gene_ontology_term_id_0",
                    "tissue_ontology_term_id_0",
                    "cell_type_ontology_term_id_2",
                ): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type
                ),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type
                ),
            },
        }

        self.assertEqual(expected, result.to_dict())


class QueryPrimaryFilterDimensionsTest(unittest.TestCase):
    def test__single_dimension__returns_all_dimension_and_terms(self):
        dim_size = 3
        with create_temp_cube(dim_size=dim_size) as cube:
            result = WmgQuery(cube).list_primary_filter_dimension_term_ids('gene_ontology_term_id')
            self.assertEquals(["gene_ontology_term_id_0",
                               "gene_ontology_term_id_1",
                               "gene_ontology_term_id_2"],
                              result)

    def test__multiple_dimensions__returns_all_dimensions_and_terms_as_tuples(self):
        dim_size = 3

        # we want disjoint set of genes across organisms, to mimic reality (each organism has its own set of genes);
        # without this filtering function, the cube would have the cross-product of organisms * genes
        def exclude(logical_coord: Tuple) -> bool:
            return (logical_coord[0], logical_coord[2]) not in \
                   {
                       ('gene_ontology_term_id_0', 'organism_ontology_term_id_0'),
                       ('gene_ontology_term_id_1', 'organism_ontology_term_id_0'),
                       ('gene_ontology_term_id_2', 'organism_ontology_term_id_1'),
                   }

        with create_temp_cube(dim_size=dim_size, exclude_logical_coord_fn=exclude) as cube:
            result = WmgQuery(cube).list_grouped_primary_filter_dimensions_term_ids(
                    'gene_ontology_term_id', 'organism_ontology_term_id')
            self.maxDiff = None
            self.assertEquals(
                    {
                        'organism_ontology_term_id_0': ['gene_ontology_term_id_0', 'gene_ontology_term_id_1'],
                        'organism_ontology_term_id_1': ['gene_ontology_term_id_2']
                    },
                    result)
