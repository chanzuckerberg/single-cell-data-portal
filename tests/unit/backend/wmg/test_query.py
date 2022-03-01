import unittest

from backend.wmg.data.query import WmgQueryCriteria, WmgQuery
from backend.wmg.data.schema import cube_non_indexed_dims
from unit.backend.wmg.fixtures.cube import create_temp_cube, all_ones_attr_values


class QueryTest(unittest.TestCase):
    def test__query_required_indexed_dims_single_value__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"]
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = WmgQuery(cube).execute(criteria)

        # sanity check the expected values of the `n` stat for each data viz point; if this fails, the cube test
        # fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture has
        # changed
        expected_cell_count_per_cell_type = dim_size ** len(
            set(cube_non_indexed_dims).difference({"cell_type_ontology_term_id"})
        )
        assert expected_cell_count_per_cell_type == 729

        expected = {
            "n_cells": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "nnz": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type),
            },
        }

        self.assertEqual(expected, result.to_dict())

    def test__query_all_indexed_dims_single_value__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_1",
            tissue_ontology_term_ids=["tissue_ontology_term_id_2"]
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = WmgQuery(cube).execute(criteria)

        # sanity check the expected values of the `n` stat for each data viz point; if this fails, the cube test
        # fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture has
        # changed
        expected_cell_count_per_cell_type = dim_size ** len(
            set(cube_non_indexed_dims).difference({"cell_type_ontology_term_id"})
        )
        assert expected_cell_count_per_cell_type == 729

        expected = {
            "n_cells": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "nnz": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type),
            },
        }

        self.assertEqual(expected, result.to_dict())

    def test__query_all_indexed_dims_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0", "gene_ontology_term_id_2"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"]
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = WmgQuery(cube).execute(criteria)

        # sanity check the expected values of the `n` stat for each data viz point; if this fails, the cube test
        # fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture has
        # changed
        expected_cell_count_per_cell_type = dim_size ** len(
            set(cube_non_indexed_dims).difference({"cell_type_ontology_term_id"})
        )
        assert expected_cell_count_per_cell_type == 729

        expected = {
            "n_cells": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "nnz": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_0"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_1"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_2"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_0"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_1"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_2"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_0"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_1"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_1",
                 "cell_type_ontology_term_id_2"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_0"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_1"): float(expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_2",
                 "cell_type_ontology_term_id_2"): float(expected_cell_count_per_cell_type),
            },
        }

        self.assertEqual(expected, result.to_dict())

    def test__query_non_indexed_dim_single_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            dataset_ids=['dataset_id_1']  # <-- non-indexed dim, single-valued
        )

        dim_size = 2
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = WmgQuery(cube).execute(criteria)

        # sanity check the expected values of the `n` stat for each data viz point; if this fails, the cube test
        # fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture has
        # changed
        expected_cell_count_per_cell_type = dim_size ** len(
            set(cube_non_indexed_dims).difference({"cell_type_ontology_term_id", "dataset_id"})
        )
        assert expected_cell_count_per_cell_type == 32

        expected = {
            "n_cells": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
            },
            "nnz": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type),
            },
        }

        self.assertEqual(expected, result.to_dict())

    def test__query_non_indexed_dim_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            dataset_ids=['dataset_id_1', 'dataset_id_0']  # <-- non-indexed dim, multi-valued
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = WmgQuery(cube).execute(criteria)

        # sanity check the expected values of the `n` stat for each data viz point; if this fails, the cube test
        # fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture has
        # changed
        expected_cell_count_per_cell_type = dim_size ** len(
            set(cube_non_indexed_dims).difference({"cell_type_ontology_term_id", "dataset_id"})
        ) * 2
        assert expected_cell_count_per_cell_type == 486

        expected = {
            "n_cells": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "nnz": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                    expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                    expected_cell_count_per_cell_type),
            },
        }

        self.assertEqual(expected, result.to_dict())

    def test__query_non_indexed_dim_single_and_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
                organism_ontology_term_id="organism_ontology_term_id_0",
                tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                assay_ontology_term_ids=['assay_ontology_term_id_1'],  # <-- non-indexed dim, single-valued
                dataset_ids=['dataset_id_1', 'dataset_id_0']  # <-- non-indexed dim, multi-valued
        )
    
        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = WmgQuery(cube).execute(criteria)
    
        # sanity check the expected values of the `n` stat for each data viz point; if this fails, the cube test
        # fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture has
        # changed
        expected_cell_count_per_cell_type = dim_size ** len(
                set(cube_non_indexed_dims).difference({"cell_type_ontology_term_id", "assay_ontology_term_id", "dataset_id"})
        ) * 1 * 2
        assert expected_cell_count_per_cell_type == 162
    
        expected = {
            "n_cells": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "nnz": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_0"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_1"): expected_cell_count_per_cell_type,
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0",
                 "cell_type_ontology_term_id_2"): expected_cell_count_per_cell_type,
            },
            "sum": {
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                        expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                        expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_0", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                        expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                        expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                        expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_1", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                        expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_0"): float(
                        expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_1"): float(
                        expected_cell_count_per_cell_type),
                ("gene_ontology_term_id_2", "tissue_ontology_term_id_0", "cell_type_ontology_term_id_2"): float(
                        expected_cell_count_per_cell_type),
            },
        }
    
        self.assertEqual(expected, result.to_dict())
