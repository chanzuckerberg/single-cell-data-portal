import unittest

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from backend.wmg.api.v1 import add_missing_combinations_to_dot_plot_matrix, rollup
from backend.wmg.data.rollup import are_cell_types_colinear, rollup_across_cell_type_descendants
from backend.wmg.data.schemas.cube_schema import expression_summary_logical_attrs


class RollupExpressionsAcrossCellTypesTest(unittest.TestCase):
    # test that the rollup function works as expected
    def test__expression_rollup_across_cell_type_descendants(self):
        # second cell type is descendant of first cell type
        # fourth cell type is descendant of third cell type
        cell_types = ["CL:0000786", "CL:0000986", "CL:0000980", "CL:0001202"]
        df = pd.DataFrame()
        df["cell_type_ontology_term_id"] = cell_types
        np.random.seed(0)
        exprs = np.random.rand(len(cell_types), 10)
        for i in range(exprs.shape[1]):
            df[i] = exprs[:, i]

        df_rollup = rollup_across_cell_type_descendants(df)
        df_expected = df.copy()
        expected_exprs = exprs.copy()
        expected_exprs[0] = exprs[0] + exprs[1]
        expected_exprs[2] = exprs[2] + exprs[3]

        for i in range(expected_exprs.shape[1]):
            df_expected[i] = expected_exprs[:, i]
        assert np.all(df_expected == df_rollup)

    def test__expression_rollup_no_descendants_overlap(self):
        # these cell types are not descendants of each other
        cell_types = ["CL:0000786", "CL:0000980"]
        df = pd.DataFrame()
        df["cell_type_ontology_term_id"] = cell_types
        np.random.seed(0)
        exprs = np.random.rand(len(cell_types), 10)
        for i in range(exprs.shape[1]):
            df[i] = exprs[:, i]

        df_rollup = rollup_across_cell_type_descendants(df)
        assert np.all(df == df_rollup)

    def test__cell_types_in_same_lineage_are_colinear(self):
        # first and second pairs of cell types are in the same lineage
        # third pair of cell types are not in the same lineage
        cell_type_pairs = [["CL:0000786", "CL:0000986"], ["CL:0000980", "CL:0001202"], ["CL:0000786", "CL:0000980"]]
        expected_colinearity = [True, True, False]
        for cell_types, expected in zip(cell_type_pairs, expected_colinearity):
            a, b = cell_types
            assert are_cell_types_colinear(a, b) == expected

    def test__missing_combinations_added_to_expression_dataframe_before_rollup(self):
        """
        test that the `add_missing_combinations_to_dot_plot_matrix` function works as expected
        this function is used to add missing (tissue, cell type) combinations to the expression dataframe
        so that the expression dataframe can be rolled up correctly.

        the cell counts dataframe has 100 cell types, 60 of which have expression data
        the goal of this test is to make sure that `add_missing_combinations_to_dot_plot_matrix`
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
        actual_dot_plot_matrix_df = add_missing_combinations_to_dot_plot_matrix(
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

        # now test the entire rollup function end-to-end, which includes `add_missing_combinations_to_dot_plot_matrix`
        dot_plot_matrix_df, cell_counts_cell_type_agg = rollup(dot_plot_matrix_df, cell_counts_cell_type_agg)

        # assert that the n_cells_tissue column is still 100
        assert (dot_plot_matrix_df["n_cells_tissue"] == 100).all()

        # assert that the added rows that could not be rescued via rollup are properly filtered out
        # these rows will still have 0 for all numeric columns, excluding n_cells_tissue
        assert (dot_plot_matrix_df["sum"] > 0).all()

    def test__gene_expression_rollup(self):
        """
        Test that the `rollup` function sums up relevant numeric values for gene expression for
        each cell type from its descendant cell types for every (tissue, gene) combination
        in the input expression dataframe.

        The input:

        1. A cell type ontology subgraph consisting of 4 cell types.
        2. 2 tissues.
        3. 2 genes.
        4. A gene expression dataframe consisting of expression numeric values for each
           (gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id) tuple.
        5. The gene expression datafram also holds total cell counts per tissue_ontology_term_id.
           This value is held in a column called `n_cells_tissue`.
        6. A cell counts dataframe that consists of cell counts for each
           (tissue_ontology_term_id, cell_type_ontology_term_id) tuple.
        7. We set the `n_cells_tissue` column values to 10_000_000 for all rows
           in the gene expression dataframe.
        8. We set all other numeric column values to `1` in the gene expression dataframe.

        The expected output:

        1. A rolled up gene expression dataframe.
        2. A rolled up cell counts dataframe.
        3. Assert that `n_cells_tissue` column value in the rolled up gene expression dataframe
           does not change for all rows because it should not be rolled up the cell type ontology
           ancestor path.
        4. Assert that other numeric column values (i.e the columns that are not `n_cells_tissue`)
           in the rolled up gene expression dataframe sum up to the correct value for each cell type
           for every (tissue, gene) combination.
        5. Assert that the cell counts in the rolled up cell counts dataframe sum up to the correct
           value for each cell type for every (tissue, cell_type) combination.
        """
        # Arrange

        # Cell Type Ontology Subgraph
        #    CL:0000127
        #    ├── CL:0000644
        #    ├── CL:0002605
        #    └── CL:0002627
        cell_counts_col_names = ["tissue_ontology_term_id", "cell_type_ontology_term_id", "n_cells_cell_type"]

        input_cell_counts_rows = [
            ["UBERON:0000955", "CL:0000127", 100000],
            ["UBERON:0000955", "CL:0000644", 8034],
            ["UBERON:0000955", "CL:0002605", 70009],
            ["UBERON:0000955", "CL:0002627", 9871],
            ["UBERON:0002113", "CL:0000127", 100000],
            ["UBERON:0002113", "CL:0000644", 8034],
            ["UBERON:0002113", "CL:0002605", 70009],
            ["UBERON:0002113", "CL:0002627", 9871],
        ]
        cell_counts_df = pd.DataFrame(input_cell_counts_rows, columns=cell_counts_col_names)
        cell_counts_df = cell_counts_df.set_index(["tissue_ontology_term_id", "cell_type_ontology_term_id"])

        gene_expr_col_names = [
            "gene_ontology_term_id",
            "tissue_ontology_term_id",
            "cell_type_ontology_term_id",
            "nnz",
            "sum",
            "n_cells_cell_type",
            "n_cells_tissue",
        ]

        input_gene_expr_rows = [
            ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 8034, 10000000],
            ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 70009, 10000000],
            ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 9871, 10000000],
            ["ENSG00000085265", "UBERON:0002113", "CL:0000644", 1, 1, 8034, 10000000],
            ["ENSG00000085265", "UBERON:0002113", "CL:0002605", 1, 1, 70009, 10000000],
            ["ENSG00000169429", "UBERON:0002113", "CL:0002627", 1, 1, 9871, 10000000],
        ]

        gene_expr_df = pd.DataFrame(input_gene_expr_rows, columns=gene_expr_col_names)

        # Act
        rolled_up_gene_expr_df, rolled_up_cell_counts_df = rollup(gene_expr_df.copy(), cell_counts_df.copy())

        # Assert
        expected_cell_counts_rows = [
            ["UBERON:0000955", "CL:0000127", 187914],
            ["UBERON:0000955", "CL:0000644", 8034],
            ["UBERON:0000955", "CL:0002605", 70009],
            ["UBERON:0000955", "CL:0002627", 9871],
            ["UBERON:0002113", "CL:0000127", 187914],
            ["UBERON:0002113", "CL:0000644", 8034],
            ["UBERON:0002113", "CL:0002605", 70009],
            ["UBERON:0002113", "CL:0002627", 9871],
        ]

        rolled_up_cell_counts_df.reset_index(inplace=True)
        rolled_up_cell_counts_df_list = rolled_up_cell_counts_df.values.tolist()

        assert rolled_up_cell_counts_df_list == expected_cell_counts_rows

        expected_gene_expr_rows = [
            ["ENSG00000085265", "UBERON:0000955", "CL:0000127", 2, 2, 78043, 10000000],
            ["ENSG00000169429", "UBERON:0000955", "CL:0000127", 1, 1, 9871, 10000000],
            ["ENSG00000085265", "UBERON:0000955", "CL:0000644", 1, 1, 8034, 10000000],
            ["ENSG00000085265", "UBERON:0000955", "CL:0002605", 1, 1, 70009, 10000000],
            ["ENSG00000169429", "UBERON:0000955", "CL:0002627", 1, 1, 9871, 10000000],
            ["ENSG00000085265", "UBERON:0002113", "CL:0000127", 2, 2, 78043, 10000000],
            ["ENSG00000169429", "UBERON:0002113", "CL:0000127", 1, 1, 9871, 10000000],
            ["ENSG00000085265", "UBERON:0002113", "CL:0000644", 1, 1, 8034, 10000000],
            ["ENSG00000085265", "UBERON:0002113", "CL:0002605", 1, 1, 70009, 10000000],
            ["ENSG00000169429", "UBERON:0002113", "CL:0002627", 1, 1, 9871, 10000000],
        ]

        rolled_up_gene_expr_df.sort_values(
            ["tissue_ontology_term_id", "cell_type_ontology_term_id", "gene_ontology_term_id"], inplace=True
        )
        rolled_up_gene_expr_list = rolled_up_gene_expr_df.values.tolist()

        assert rolled_up_gene_expr_list == expected_gene_expr_rows
