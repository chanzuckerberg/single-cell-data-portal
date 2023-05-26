"""This module tests the rollup functions called directly and indirectly by the WMG API.

In detail, this module tests the public and private functions defined in `backend.wmg.api.common.rollup` module.
"""
import unittest

import pandas as pd
from pandas import DataFrame
from pandas.testing import assert_frame_equal
from parameterized import parameterized

from backend.wmg.api.common.rollup import _add_missing_combinations_to_gene_expression_df_for_rollup, rollup
from backend.wmg.data.schemas.cube_schema import expression_summary_logical_attrs


def _create_cell_counts_df_helper(cell_counts_rows: list[list], columns: list[str], index_cols: list[str]) -> DataFrame:
    cell_counts_df = pd.DataFrame(cell_counts_rows, columns=columns)
    cell_counts_df = cell_counts_df.set_index(index_cols, verify_integrity=True)
    return cell_counts_df


def _create_gene_expression_df_helper(gene_expr_rows: list[list], columns: list[str]) -> DataFrame:
    gene_expr_df = pd.DataFrame(gene_expr_rows, columns=columns)
    return gene_expr_df


def _cell_counts_df_without_compare_dim(cell_counts_rows: list[list]) -> DataFrame:
    cell_counts_col_names = ["tissue_ontology_term_id", "cell_type_ontology_term_id", "n_cells_cell_type"]
    cell_counts_index_col_names = ["tissue_ontology_term_id", "cell_type_ontology_term_id"]
    return _create_cell_counts_df_helper(
        cell_counts_rows, columns=cell_counts_col_names, index_cols=cell_counts_index_col_names
    )


def _cell_counts_df_with_ethnicity_compare_dim(cell_counts_rows: list[list]) -> DataFrame:
    cell_counts_col_names = [
        "tissue_ontology_term_id",
        "cell_type_ontology_term_id",
        "self_reported_ethnicity_ontology_term_id",
        "n_cells_cell_type",
    ]

    cell_counts_index_col_names = [
        "tissue_ontology_term_id",
        "cell_type_ontology_term_id",
        "self_reported_ethnicity_ontology_term_id",
    ]

    return _create_cell_counts_df_helper(
        cell_counts_rows, columns=cell_counts_col_names, index_cols=cell_counts_index_col_names
    )


def _gene_expression_df_without_compare_dim(gene_expr_rows: list[list]) -> DataFrame:
    gene_expr_col_names = [
        "gene_ontology_term_id",
        "tissue_ontology_term_id",
        "cell_type_ontology_term_id",
        "nnz",
        "sum",
        "n_cells_cell_type",
        "n_cells_tissue",
    ]
    return _create_gene_expression_df_helper(gene_expr_rows, columns=gene_expr_col_names)


def _gene_expression_df_with_ethnicity_compare_dim(gene_expr_rows: list[list]) -> DataFrame:
    gene_expr_col_names = [
        "gene_ontology_term_id",
        "tissue_ontology_term_id",
        "cell_type_ontology_term_id",
        "self_reported_ethnicity_ontology_term_id",
        "nnz",
        "sum",
        "n_cells_cell_type",
        "n_cells_tissue",
    ]

    return _create_gene_expression_df_helper(gene_expr_rows, columns=gene_expr_col_names)


class TestHighLevelRollupFunction(unittest.TestCase):
    """
    Test that the `rollup` function correctly accumulates (or rolls up) gene-expression
    values FOR EACH expressed gene and cell count values up the cell type ANCESTOR paths
    grouped by (tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>).

    The input:

    1. A cell type ontology subgraph consisting of 4 cell types.

    CL:0000127
    ├── CL:0000644
    ├── CL:0002605
    └── CL:0002627

    2. A gene expression dataframe consisting of expression numeric values for each
    (gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>) tuple.
    3. The gene expression dataframe also holds total cell counts per tissue_ontology_term_id.
    This value is held in a column called `n_cells_tissue`.
    6. A cell counts dataframe that consists of cell counts for each
    (tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>) tuple.
    7. We set the `n_cells_tissue` column values to 1000 for all rows
    in the gene expression dataframe.
    8. We set all other numeric column values to `1` in the gene expression dataframe.

    The expected output:

    1. A rolled up gene expression dataframe.
    2. A rolled up cell counts dataframe.
    3. Assert that `n_cells_tissue` column value in the rolled up gene expression dataframe
    does not change for all rows because it should not be rolled up the cell type ontology
    ancestor paths.
    4. Assert that other numeric column values (i.e the columns that are not `n_cells_tissue`)
    in the rolled up gene expression dataframe hold the correct rolled up values.
    5. Assert that the cell counts in the rolled up cell counts dataframe hold the correct
    rolled up values.
    """

    @staticmethod
    def _rollup_testcases():
        """
        Testcases for the `rollup` function.

        An important note about how the expected values are laid out in the testcases:

        1. Expected values for rows in the rolled up cell counts dataframe are sorted by
        (tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>)

        2. Expected values for rows in the rolled up gene expression dataframe are sorted by
        (tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>, gene_ontology_term_id)
        """
        test_1 = {
            "name": "no_compare_dim_all_tissues_have_all_cell_types",
            "input_cell_counts": [
                ["UBERON:0000955", "CL:0000127", 300],
                ["UBERON:0000955", "CL:0000644", 70],
                ["UBERON:0000955", "CL:0002605", 80],
                ["UBERON:0000955", "CL:0002627", 90],
                ["UBERON:0002113", "CL:0000127", 300],
                ["UBERON:0002113", "CL:0000644", 70],
                ["UBERON:0002113", "CL:0002605", 80],
                ["UBERON:0002113", "CL:0002627", 90],
            ],
            "expected_rolled_up_cell_counts": [
                ["UBERON:0000955", "CL:0000127", 540],
                ["UBERON:0000955", "CL:0000644", 70],
                ["UBERON:0000955", "CL:0002605", 80],
                ["UBERON:0000955", "CL:0002627", 90],
                ["UBERON:0002113", "CL:0000127", 540],
                ["UBERON:0002113", "CL:0000644", 70],
                ["UBERON:0002113", "CL:0002605", 80],
                ["UBERON:0002113", "CL:0002627", 90],
            ],
            "input_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0002113", "CL:0002627", 1, 1, 90, 1000],
            ],
            "expected_rolled_up_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000127", 2, 2, 150, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0000127", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0000127", 2, 2, 150, 1000],
                ["ENSG00000169429", "UBERON:0002113", "CL:0000127", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0002113", "CL:0002627", 1, 1, 90, 1000],
            ],
        }

        test_2 = {
            "name": "no_compare_dim_one_ancestor_cell_type_missing_in_one_tissue_but_exists_in_all_others",
            # Tissue: "UBERON:0000955" MISSING cell type: "CL:0000127" in input cell counts
            "input_cell_counts": [
                ["UBERON:0000955", "CL:0000644", 70],
                ["UBERON:0000955", "CL:0002605", 80],
                ["UBERON:0000955", "CL:0002627", 90],
                ["UBERON:0002113", "CL:0000127", 300],
                ["UBERON:0002113", "CL:0000644", 70],
                ["UBERON:0002113", "CL:0002605", 80],
                ["UBERON:0002113", "CL:0002627", 90],
            ],
            # cell count for cell type: "CL:0000127" in Tissue: "UBERON:0000955" GETS AGGREGATED
            # count because Tissue "UBERON:0000955" CONTAINS cell counts for
            # descendants of cell type: "CL:0000127"
            "expected_rolled_up_cell_counts": [
                ["UBERON:0000955", "CL:0000127", 240],
                ["UBERON:0000955", "CL:0000644", 70],
                ["UBERON:0000955", "CL:0002605", 80],
                ["UBERON:0000955", "CL:0002627", 90],
                ["UBERON:0002113", "CL:0000127", 540],
                ["UBERON:0002113", "CL:0000644", 70],
                ["UBERON:0002113", "CL:0002605", 80],
                ["UBERON:0002113", "CL:0002627", 90],
            ],
            "input_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0002113", "CL:0002627", 1, 1, 90, 1000],
            ],
            "expected_rolled_up_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000127", 2, 2, 150, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0000127", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0000127", 2, 2, 150, 1000],
                ["ENSG00000169429", "UBERON:0002113", "CL:0000127", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0002113", "CL:0002627", 1, 1, 90, 1000],
            ],
        }

        test_3 = {
            "name": "no_compare_dim_gene_expressed_in_one_tissue_but_not_other",
            "input_cell_counts": [
                ["UBERON:0000955", "CL:0000127", 300],
                ["UBERON:0000955", "CL:0000644", 70],
                ["UBERON:0000955", "CL:0002605", 80],
                ["UBERON:0000955", "CL:0002627", 90],
                ["UBERON:0002113", "CL:0000127", 300],
                ["UBERON:0002113", "CL:0000644", 70],
                ["UBERON:0002113", "CL:0002605", 80],
                ["UBERON:0002113", "CL:0002627", 90],
            ],
            "expected_rolled_up_cell_counts": [
                ["UBERON:0000955", "CL:0000127", 540],
                ["UBERON:0000955", "CL:0000644", 70],
                ["UBERON:0000955", "CL:0002605", 80],
                ["UBERON:0000955", "CL:0002627", 90],
                ["UBERON:0002113", "CL:0000127", 540],
                ["UBERON:0002113", "CL:0000644", 70],
                ["UBERON:0002113", "CL:0002605", 80],
                ["UBERON:0002113", "CL:0002627", 90],
            ],
            # Gene "ENSG00000169429" expressed in Tissue "UBERON:0000955" but not expressed
            # in Tissue "UBERON:0002113"
            "input_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0002605", 1, 1, 80, 1000],
            ],
            "expected_rolled_up_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000127", 2, 2, 150, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0000127", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0000127", 2, 2, 150, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0002113", "CL:0002605", 1, 1, 80, 1000],
            ],
        }

        test_4 = {
            "name": "no_compare_dim_one_of_the_tissues_has_no_gene_expressions_at_all",
            "input_cell_counts": [
                ["UBERON:0000955", "CL:0000127", 300],
                ["UBERON:0000955", "CL:0000644", 70],
                ["UBERON:0000955", "CL:0002605", 80],
                ["UBERON:0000955", "CL:0002627", 90],
                ["UBERON:0002113", "CL:0000127", 300],
                ["UBERON:0002113", "CL:0000644", 70],
                ["UBERON:0002113", "CL:0002605", 80],
                ["UBERON:0002113", "CL:0002627", 90],
            ],
            "expected_rolled_up_cell_counts": [
                ["UBERON:0000955", "CL:0000127", 540],
                ["UBERON:0000955", "CL:0000644", 70],
                ["UBERON:0000955", "CL:0002605", 80],
                ["UBERON:0000955", "CL:0002627", 90],
                ["UBERON:0002113", "CL:0000127", 540],
                ["UBERON:0002113", "CL:0000644", 70],
                ["UBERON:0002113", "CL:0002605", 80],
                ["UBERON:0002113", "CL:0002627", 90],
            ],
            # Tissue issue "UBERON:0002113" has no gene expressions
            "input_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 90, 1000],
            ],
            "expected_rolled_up_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000127", 2, 2, 150, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0000127", 1, 1, 90, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 80, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 90, 1000],
            ],
        }

        test_5 = {
            "name": "with_ethnicity_compare_dim_on_single_tissue",
            "input_cell_counts": [
                ["UBERON:0000955", "CL:0000127", "unknown", 300],
                ["UBERON:0000955", "CL:0000644", "unknown", 70],
                ["UBERON:0000955", "CL:0002605", "HANCESTRO:0005", 10],
                ["UBERON:0000955", "CL:0002605", "HANCESTRO:0008", 30],
                ["UBERON:0000955", "CL:0002605", "multiethnic", 40],
                ["UBERON:0000955", "CL:0002627", "HANCESTRO:0006", 10],
                ["UBERON:0000955", "CL:0002627", "HANCESTRO:0008", 20],
                ["UBERON:0000955", "CL:0002627", "multiethnic", 30],
                ["UBERON:0000955", "CL:0002627", "unknown", 40],
            ],
            "expected_rolled_up_cell_counts": [
                ["UBERON:0000955", "CL:0000127", "HANCESTRO:0005", 10],
                ["UBERON:0000955", "CL:0000127", "HANCESTRO:0006", 10],
                ["UBERON:0000955", "CL:0000127", "HANCESTRO:0008", 50],
                ["UBERON:0000955", "CL:0000127", "multiethnic", 70],
                ["UBERON:0000955", "CL:0000127", "unknown", 410],
                ["UBERON:0000955", "CL:0000644", "unknown", 70],
                ["UBERON:0000955", "CL:0002605", "HANCESTRO:0005", 10],
                ["UBERON:0000955", "CL:0002605", "HANCESTRO:0008", 30],
                ["UBERON:0000955", "CL:0002605", "multiethnic", 40],
                ["UBERON:0000955", "CL:0002627", "HANCESTRO:0006", 10],
                ["UBERON:0000955", "CL:0002627", "HANCESTRO:0008", 20],
                ["UBERON:0000955", "CL:0002627", "multiethnic", 30],
                ["UBERON:0000955", "CL:0002627", "unknown", 40],
            ],
            "input_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", "unknown", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", "HANCESTRO:0008", 1, 1, 30, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", "multiethnic", 1, 1, 40, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002627", "HANCESTRO:0008", 1, 1, 20, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", "multiethnic", 1, 1, 30, 1000],
            ],
            "expected_rolled_up_gene_expression": [
                ["ENSG00000085265", "UBERON:0000955", "CL:0000127", "HANCESTRO:0008", 2, 2, 50, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0000127", "multiethnic", 1, 1, 40, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0000127", "multiethnic", 1, 1, 30, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0000127", "unknown", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0000644", "unknown", 1, 1, 70, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", "HANCESTRO:0008", 1, 1, 30, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002605", "multiethnic", 1, 1, 40, 1000],
                ["ENSG00000085265", "UBERON:0000955", "CL:0002627", "HANCESTRO:0008", 1, 1, 20, 1000],
                ["ENSG00000169429", "UBERON:0000955", "CL:0002627", "multiethnic", 1, 1, 30, 1000],
            ],
        }

        return [
            (
                test_1["name"],
                _cell_counts_df_without_compare_dim(test_1["input_cell_counts"]),
                _cell_counts_df_without_compare_dim(test_1["expected_rolled_up_cell_counts"]),
                _gene_expression_df_without_compare_dim(test_1["input_gene_expression"]),
                _gene_expression_df_without_compare_dim(test_1["expected_rolled_up_gene_expression"]),
            ),
            (
                test_2["name"],
                _cell_counts_df_without_compare_dim(test_2["input_cell_counts"]),
                _cell_counts_df_without_compare_dim(test_2["expected_rolled_up_cell_counts"]),
                _gene_expression_df_without_compare_dim(test_2["input_gene_expression"]),
                _gene_expression_df_without_compare_dim(test_2["expected_rolled_up_gene_expression"]),
            ),
            (
                test_3["name"],
                _cell_counts_df_without_compare_dim(test_3["input_cell_counts"]),
                _cell_counts_df_without_compare_dim(test_3["expected_rolled_up_cell_counts"]),
                _gene_expression_df_without_compare_dim(test_3["input_gene_expression"]),
                _gene_expression_df_without_compare_dim(test_3["expected_rolled_up_gene_expression"]),
            ),
            (
                test_4["name"],
                _cell_counts_df_without_compare_dim(test_4["input_cell_counts"]),
                _cell_counts_df_without_compare_dim(test_4["expected_rolled_up_cell_counts"]),
                _gene_expression_df_without_compare_dim(test_4["input_gene_expression"]),
                _gene_expression_df_without_compare_dim(test_4["expected_rolled_up_gene_expression"]),
            ),
            (
                test_5["name"],
                _cell_counts_df_with_ethnicity_compare_dim(test_5["input_cell_counts"]),
                _cell_counts_df_with_ethnicity_compare_dim(test_5["expected_rolled_up_cell_counts"]),
                _gene_expression_df_with_ethnicity_compare_dim(test_5["input_gene_expression"]),
                _gene_expression_df_with_ethnicity_compare_dim(test_5["expected_rolled_up_gene_expression"]),
            ),
        ]

    @parameterized.expand(_rollup_testcases)
    def test__rollup(self, _, input_cell_counts_df, expected_cell_counts_df, input_gene_expr_df, expected_gene_expr_df):
        # Arrange
        cell_counts_df_index_list = list(input_cell_counts_df.index.names)

        # Act

        # Note that we are creating copies of the input dataframes before passing them as
        # arguments to the `rollup` function so that if the `rollup` function mutates the
        # argument values, the input to the test is not affected.
        rolled_up_gene_expr_df, rolled_up_cell_counts_df = rollup(
            input_gene_expr_df.copy(), input_cell_counts_df.copy()
        )

        # Assert
        rolled_up_cell_counts_df.reset_index(inplace=True)
        expected_cell_counts_df.reset_index(inplace=True)

        assert_frame_equal(
            rolled_up_cell_counts_df.reset_index(drop=True), expected_cell_counts_df.reset_index(drop=True)
        )

        # sort the rolled up gene expression dataframe so that the correct rows are compared with
        # the expected gene expression rows in the assert call
        sort_columns_for_rolled_gene_expr_df = list(cell_counts_df_index_list) + ["gene_ontology_term_id"]
        rolled_up_gene_expr_df.sort_values(sort_columns_for_rolled_gene_expr_df, inplace=True)

        assert_frame_equal(rolled_up_gene_expr_df.reset_index(drop=True), expected_gene_expr_df.reset_index(drop=True))


class TestHighLevelRollupHelperFunctions(unittest.TestCase):
    def test__add_missing_combinations_to_gene_expression_df_for_rollup(self):
        """
        test that the `_add_missing_combinations_to_gene_expression_df_for_rollup` function works as expected
        this function is used to add missing (tissue, cell type) combinations to the expression dataframe
        so that the expression dataframe can be rolled up correctly.

        the cell counts dataframe has 100 cell types, 60 of which have expression data
        the goal of this test is to make sure that `_add_missing_combinations_to_gene_expression_df_for_rollup`
        adds the missing (tissue, cell type) combinations for each gene to th expression dataframe.

        there will be (100 cell types) * (2 tissues) = 200 rows in the cell counts dataframe
        there will be (60 cell types) * (2 tissues) * (3 genes) = 360 rows in the expression dataframe
        there will be (100 cell types) * (2 tissues) * (3 genes) = 600 rows in the expected expression dataframe


        the input expression dataframe will have 1 for all numeric columns, excluding n_cells_tissue.
        because we don't want to roll up the n_cells_tissue column, we set it to 100 for all rows and check
        that it stays 100 after adding the missing combinations.

        the expected expression dataframe will have the same rows as the input dataframe, with new rows added
        that have 0 for all numeric columns, excluding n_cells_tissue, which stays 100.
        """
        num_cell_types = 100
        num_cell_types_with_expression = 60
        num_tissues = 2
        num_genes = 3

        cell_types_cell_counts_df = [f"cell_type_ontology_term_id_{i}" for i in range(num_cell_types)] * num_tissues
        tissue_cell_counts_df = sum([[f"tissue_ontology_term_id_{i}"] * num_cell_types for i in range(num_tissues)], [])

        cell_types_dot_plot_df = cell_types_cell_counts_df[:num_cell_types_with_expression] * num_tissues * num_genes
        tissue_dot_plot_df = (
            sum([[f"tissue_ontology_term_id_{i}"] * num_cell_types_with_expression for i in range(num_tissues)], [])
            * num_genes
        )
        genes_dot_plot_df = sum(
            [[f"gene_ontology_term_id_{i}"] * num_cell_types_with_expression * num_tissues for i in range(num_genes)],
            [],
        )

        cell_counts_cell_type_agg = pd.DataFrame()
        cell_counts_cell_type_agg["cell_type_ontology_term_id"] = cell_types_cell_counts_df
        cell_counts_cell_type_agg["tissue_ontology_term_id"] = tissue_cell_counts_df
        cell_counts_cell_type_agg["n_cells_cell_type"] = 1

        dot_plot_matrix_df = pd.DataFrame()
        dot_plot_matrix_df["cell_type_ontology_term_id"] = cell_types_dot_plot_df
        dot_plot_matrix_df["tissue_ontology_term_id"] = tissue_dot_plot_df
        dot_plot_matrix_df["gene_ontology_term_id"] = genes_dot_plot_df
        for attr in expression_summary_logical_attrs:
            dot_plot_matrix_df[attr.name] = 1
        dot_plot_matrix_df["n_cells_cell_type"] = 1
        dot_plot_matrix_df["n_cells_tissue"] = 100

        expected_added_cell_types = (
            [f"cell_type_ontology_term_id_{i}" for i in range(num_cell_types_with_expression, num_cell_types)]
            * num_tissues
            * num_genes
        )
        expected_added_tissues = (
            sum(
                [
                    [f"tissue_ontology_term_id_{i}"] * (num_cell_types - num_cell_types_with_expression)
                    for i in range(num_tissues)
                ],
                [],
            )
            * num_genes
        )
        expected_added_genes = sum(
            [
                [f"gene_ontology_term_id_{i}"] * (num_cell_types - num_cell_types_with_expression) * num_tissues
                for i in range(num_genes)
            ],
            [],
        )
        expected_added_dataframe = pd.DataFrame()
        expected_added_dataframe["cell_type_ontology_term_id"] = expected_added_cell_types
        expected_added_dataframe["tissue_ontology_term_id"] = expected_added_tissues
        expected_added_dataframe["gene_ontology_term_id"] = expected_added_genes
        for attr in expression_summary_logical_attrs:
            expected_added_dataframe[attr.name] = 0
        expected_added_dataframe["n_cells_cell_type"] = 0
        expected_added_dataframe["n_cells_tissue"] = 100

        expected_dot_plot_matrix_df = pd.concat((dot_plot_matrix_df, expected_added_dataframe), axis=0)

        # group by to get cell_counts_cell_type_agg into the right format (with multi-index)
        cell_counts_cell_type_agg = cell_counts_cell_type_agg.groupby(
            ["cell_type_ontology_term_id", "tissue_ontology_term_id"]
        ).first()
        actual_dot_plot_matrix_df = _add_missing_combinations_to_gene_expression_df_for_rollup(
            dot_plot_matrix_df, cell_counts_cell_type_agg
        )

        # sort both dataframes to ensure comparable ordering
        actual_dot_plot_matrix_df = actual_dot_plot_matrix_df.sort_values(
            ["cell_type_ontology_term_id", "tissue_ontology_term_id", "gene_ontology_term_id"]
        ).reset_index(drop=True)
        expected_dot_plot_matrix_df = expected_dot_plot_matrix_df.sort_values(
            ["cell_type_ontology_term_id", "tissue_ontology_term_id", "gene_ontology_term_id"]
        ).reset_index(drop=True)

        # assert equality
        assert_frame_equal(actual_dot_plot_matrix_df, expected_dot_plot_matrix_df)
