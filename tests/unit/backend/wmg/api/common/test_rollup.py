"""This module tests the rollup functions called directly and indirectly by the WMG API.

In detail, this module tests the public and private functions defined in `backend.wmg.api.common.rollup` module.
"""
import unittest

import pandas as pd
from pandas import DataFrame
from pandas.testing import assert_frame_equal
from parameterized import parameterized

from backend.wmg.api.common.rollup import rollup


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
            rolled_up_cell_counts_df.reset_index(drop=True),
            expected_cell_counts_df.reset_index(drop=True),
            check_dtype=False,
        )

        # sort the rolled up gene expression dataframe so that the correct rows are compared with
        # the expected gene expression rows in the assert call
        sort_columns_for_rolled_gene_expr_df = list(cell_counts_df_index_list) + ["gene_ontology_term_id"]
        rolled_up_gene_expr_df.sort_values(sort_columns_for_rolled_gene_expr_df, inplace=True)

        assert_frame_equal(
            rolled_up_gene_expr_df.reset_index(drop=True),
            expected_gene_expr_df.reset_index(drop=True),
            check_dtype=False,
        )
