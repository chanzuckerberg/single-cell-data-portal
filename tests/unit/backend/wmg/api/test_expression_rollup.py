import unittest

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from backend.wmg.api.v1 import add_missing_combinations_to_dot_plot_matrix, rollup
from backend.wmg.data.rollup import rollup_across_cell_type_descendants
from backend.wmg.data.schemas.cube_schema import expression_summary_logical_attrs


class RollupExpressionsAcrossCellTypesTest(unittest.TestCase):
    # test that the rollup function works as expected
    def test__expression_rollup(self):
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
        test_dot_plot_matrix_df = add_missing_combinations_to_dot_plot_matrix(
            dot_plot_matrix_df, cell_counts_cell_type_agg
        )

        # sort both dataframes to ensure comparable ordering
        test_dot_plot_matrix_df = test_dot_plot_matrix_df.sort_values(
            ["cell_type_ontology_term_id", "tissue_ontology_term_id", "gene_ontology_term_id"]
        ).reset_index(drop=True)
        expected_dot_plot_matrix_df = expected_dot_plot_matrix_df.sort_values(
            ["cell_type_ontology_term_id", "tissue_ontology_term_id", "gene_ontology_term_id"]
        ).reset_index(drop=True)

        # assert equality
        assert_frame_equal(test_dot_plot_matrix_df, expected_dot_plot_matrix_df)

        # now test the entire rollup function end-to-end, which includes `add_missing_combinations_to_dot_plot_matrix`
        dot_plot_matrix_df, cell_counts_cell_type_agg = rollup(dot_plot_matrix_df, cell_counts_cell_type_agg)

        # assert that the n_cells_tissue column is still 100
        assert (dot_plot_matrix_df["n_cells_tissue"] == 100).all()

        # assert that the added rows that could not be rescued via rollup are properly filtered out
        # these rows will still have 0 for all numeric columns, excluding n_cells_tissue
        assert (dot_plot_matrix_df["sum"] > 0).all()
