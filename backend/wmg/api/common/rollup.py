"""This module contains the implementation of the gene expression and cell count rollup feature.

The API public methods call the public methods in this module to perform the rollup operations.
"""

import itertools
from typing import Tuple

import numpy as np
import pandas as pd
from pandas import DataFrame

from backend.common.utils.rollup import rollup_across_cell_type_descendants, rollup_across_cell_type_descendants_array

######################### PUBLIC FUNCTIONS IN ALPHABETIC ORDER ##################################


def rollup(gene_expression_df, cell_counts_grouped_df) -> Tuple[DataFrame, DataFrame]:
    """
    Accumulates (or rolls up) cell count values and gene-expression values FOR EACH expressed gene
    up the cell type ANCESTOR paths grouped by the ontology term IDs in the multi-index of the
    input cell counts dataframe.

    Parameters
    ----------
    gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe containing the dimensions across which the numeric columns will be
        aggregated.

    cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe containing the dimensions across which the cell count values will be
        aggregated.

    Returns
    -------
    rolled_up_gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe with the same columns as the input gene expression dataframe,
        and likely greater size than the input gene expression dataframe, but with the numeric
        columns aggregated across the cell type's descendants.

    rolled_up_cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe with the same columns as the input cell counts dataframe,
        and likely greater size than the input cell counts dataframe, but with the cell count
        values aggregated across the cell type's descendants.
    """
    # An implementation detail note:

    # The order of operation is important: The cell counts are rolled up first to produce a the
    # `rolled_up_cell_counts_grouped_df` dataframe. Then this `rolled_up_cell_counts_grouped_df` is
    # an input into the operation to roll up gene expression values.

    # This is done for efficiency reasons where the `rolled_up_cell_counts_grouped_df` is made sparse
    # by the first operation that rolls up the cell counts. Having a sparse
    # `rolled_up_cell_counts_grouped_df` significantly improves the running time and memory footprint of
    # the second operation: rolling up the gene expression values.

    rolled_up_cell_counts_grouped_df = _rollup_cell_counts(cell_counts_grouped_df)
    rolled_up_gene_expression_df = _rollup_gene_expression(gene_expression_df, rolled_up_cell_counts_grouped_df)
    return rolled_up_gene_expression_df, rolled_up_cell_counts_grouped_df


######################### PRIVATE FUNCTIONS IN ALPHABETIC ORDER ##################################


def _rollup_gene_expression(gene_expression_df, universal_set_cell_counts_df) -> DataFrame:
    """
    Augments the input gene expression dataframe to include
    (gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>)
    combinations for which numeric expression values should be aggregated during the rollup operation.

    Parameters
    ----------
    gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe.

    universal_set_cell_counts_df : pandas DataFrame
        Multi-indexed cell counts dataframe that contains "The Universal Set Of GroupBy Ontology Terms".

    Returns
    -------
    gene_expression_with_missing_combos_df : pandas DataFrame
        Tidy gene expression dataframe with the same columns as the input gene expression dataframe,
        and likely greater size than the input gene expression dataframe that includes combinations
        for which numeric values should be aggregated during the rollup operation.
    """
    if gene_expression_df.shape[0] == 0:
        return gene_expression_df

    numeric_columns = list(
        gene_expression_df.columns[[np.issubdtype(dtype, np.number) for dtype in gene_expression_df.dtypes]]
    )

    group_by_terms = list(universal_set_cell_counts_df.index.names)
    available_combinations = set(universal_set_cell_counts_df.index.values)

    numeric_columns.remove("n_cells_tissue")

    pivoted = gene_expression_df.pivot_table(
        index=group_by_terms, columns=["gene_ontology_term_id"], fill_value=0, values=numeric_columns, aggfunc="sum"
    )
    genes = np.array(pivoted["sum"].columns.tolist())
    missing_combinations = list(available_combinations.difference(pivoted.index.values))
    dense_tables = np.stack(
        [
            np.vstack((pivoted[col].values, np.zeros((len(missing_combinations), len(genes)))))
            for col in numeric_columns
        ],
        axis=2,
    )

    multi_index = pd.MultiIndex.from_tuples(pivoted.index.tolist() + missing_combinations, names=group_by_terms)
    group_by_terms.remove("cell_type_ontology_term_id")
    group_by_terms = ["cell_type_ontology_term_id"] + group_by_terms

    multi_index = multi_index.reorder_levels(group_by_terms)

    multi_index_array = np.vstack(multi_index.values)

    cell_types_for_rollup = np.char.add(multi_index_array[:, 0], ";;")
    for col in range(1, multi_index_array.shape[1]):
        if col > 1:
            cell_types_for_rollup = np.char.add(cell_types_for_rollup, "--")
        cell_types_for_rollup = np.char.add(cell_types_for_rollup, multi_index_array[:, col])

    dense_tables = rollup_across_cell_type_descendants_array(dense_tables, cell_types_for_rollup)

    row, col = dense_tables[:, :, 0].nonzero()

    new_df = pd.DataFrame(index=multi_index[row], data=genes[col], columns=["gene_ontology_term_id"]).reset_index()
    for i, col_name in enumerate(numeric_columns):
        new_df[col_name] = dense_tables[row, col, i]

    # get the total number of cells per tissue to populate the n_cells_tissue in the added entries
    n_cells_tissue = gene_expression_df.groupby("tissue_ontology_term_id").first()["n_cells_tissue"]

    new_df["n_cells_tissue"] = n_cells_tissue[new_df["tissue_ontology_term_id"]].values
    return new_df[gene_expression_df.columns]


def _build_cell_count_groups_universal_set(cell_counts_grouped_df) -> DataFrame:
    """
    Constructs a dataframe that contains all valid combination of
    (tissue, cell_type, <compare_dimension>) for which aggregation should be performed.
    We call this set of all valid combinations "The Universal Set Of GroupBy Ontology Terms".

    The combinations are encoded as a multi-index in the input cell counts dataframe.
    The "compare dimension" is optional and may not be present in the multi-index; In that case,
    the combinations considered are (tissue, cell_type) pairs.

    In this implementation, "The Universal Set Of GroupBy Ontology Terms" is constructed by
    computing the cartesian product of the (tissue, cell_type, <compare_dimension>) values
    in the input cell counts dataframe.

    Parameters
    ----------
    cell_counts_grouped_df : pandas DataFrame
        MultiIndexed cell counts dataframe containing the dimensions across which the cell count values exist.

    Returns
    -------
    universal_set_cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe that contains "The Universal Set Of GroupBy Ontology Terms"
        with the same columns as the input cell counts dataframe, and likely greater size than the
        input cell counts dataframe.
    """
    cartesian_product = list(
        itertools.product(
            *[cell_counts_grouped_df.index.levels[i] for i in range(len(cell_counts_grouped_df.index.names))]
        )
    )
    cartesian_product_index = pd.Index(cartesian_product)
    cartesian_product_index.set_names(cell_counts_grouped_df.index.names, inplace=True)
    universal_set_cell_counts_grouped_df = pd.DataFrame(index=cartesian_product_index)
    for c in cell_counts_grouped_df.columns:
        universal_set_cell_counts_grouped_df[c] = 0
        universal_set_cell_counts_grouped_df[c][cell_counts_grouped_df.index] = cell_counts_grouped_df[c]
    return universal_set_cell_counts_grouped_df


def _rollup_cell_counts(cell_counts_grouped_df) -> DataFrame:
    """
    Roll up cell count values across cell type descendants in the input cell counts dataframe.

    Accumulate cell count values up the cell type ontology ancestor paths for every
    (tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>) combination in the
    input cell counts dataframe. The compare dimension is an optional
    field and in the case that it is missing, then accumulation happens for every
    (tissue_ontology_term_id, cell_type_ontology_term_id) combination.

    For example:

    If (T1, C1) has cell counts in the input cell counts dataframe, and the tissue labeled T1
    has another cell labeled C2 where C2 is an ANCESTOR of C1, then the rolled up cell counts dataframe
    must include CUMULATIVE cell counts values for the combination (T1, C2) by adding in the
    cell count for (T1, C1).

    Parameters
    ----------
    cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe containing the dimensions across which the cell count values will be
        aggregated.

    Returns
    -------
    rolled_up_cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe with the same columns as the input cell counts dataframe,
        and likely greater size than the input cell counts dataframe, but with the cell count
        values aggregated across the cell type's descendants.
    """
    rolled_up_cell_counts_grouped_df = cell_counts_grouped_df

    if cell_counts_grouped_df.shape[0] > 0:
        universal_set_cell_counts_grouped_df = _build_cell_count_groups_universal_set(cell_counts_grouped_df)
        index_names = universal_set_cell_counts_grouped_df.index.names

        universal_set_cell_counts_grouped_df = universal_set_cell_counts_grouped_df.reset_index()

        # rollup cell counts across cell type descendants
        rolled_up_cell_counts_grouped_df = rollup_across_cell_type_descendants(universal_set_cell_counts_grouped_df)

        rolled_up_cell_counts_grouped_df = rolled_up_cell_counts_grouped_df[
            rolled_up_cell_counts_grouped_df["n_cells_cell_type"] > 0
        ]
        rolled_up_cell_counts_grouped_df.set_index(index_names, inplace=True)

    return rolled_up_cell_counts_grouped_df
