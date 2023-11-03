"""This module contains the implementation of the gene expression and cell count rollup feature.

The API public methods call the public methods in this module to perform the rollup operations.
"""

import itertools
from typing import Tuple

import numpy as np
import pandas as pd
from ddtrace import tracer
from pandas import DataFrame

from backend.common.utils.rollup import (
    are_cell_types_not_redundant_nodes,
    rollup_across_cell_type_descendants,
    rollup_across_cell_type_descendants_array,
)

######################### PUBLIC FUNCTIONS IN ALPHABETIC ORDER ##################################


@tracer.wrap(name="rollup", service="wmg-api", resource="query", span_type="wmg-api")
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


@tracer.wrap(name="_rollup_gene_expression", service="wmg-api", resource="rollup", span_type="wmg-api")
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
    # if the gene expression dataframe is empty, then there is nothing to roll up
    if gene_expression_df.shape[0] == 0:
        return gene_expression_df

    # get the numeric columns in the gene expression dataframe
    numeric_columns = list(
        gene_expression_df.columns[[np.issubdtype(dtype, np.number) for dtype in gene_expression_df.dtypes]]
    )

    # get the available combinations of (tissue, cell_type, <compare_dimension>) in the input cell counts dataframe
    group_by_terms = list(universal_set_cell_counts_df.index.names)
    available_combinations = set(universal_set_cell_counts_df.index.values)

    # we do not wish to roll up the n_cells_tissue column
    numeric_columns.remove("n_cells_tissue")

    # pivot the gene expression dataframe to get a dense 2D array of numeric values (sum, nnz, sqsum)
    # rows are the (tissue, cell_type, <compare_dimension>) combinations and columns are the genes
    pivoted = gene_expression_df.pivot_table(
        index=group_by_terms, columns=["gene_ontology_term_id"], fill_value=0, values=numeric_columns, aggfunc="sum"
    )
    genes = np.array(pivoted["sum"].columns.tolist())
    missing_combinations = list(available_combinations.difference(pivoted.index.values))

    # vertically stack empty arrays corresponding to the missing combinations, and then
    # stack each 2D array corresponding to sum, nnz, or sqsum into a 3D array
    # this 3D array will be rolled up along the first dimension
    dense_tables = np.stack(
        [
            np.vstack((pivoted[col].values, np.zeros((len(missing_combinations), len(genes)))))
            for col in numeric_columns
        ],
        axis=2,
    )
    # get the multi-index (rows) for the pivoted tables and add the missing combinations
    multi_index = pd.MultiIndex.from_tuples(pivoted.index.tolist() + missing_combinations, names=group_by_terms)

    # reorder the multi-index to have cell_type_ontology_term_id as the first level
    group_by_terms.remove("cell_type_ontology_term_id")
    group_by_terms = ["cell_type_ontology_term_id"] + group_by_terms
    multi_index = multi_index.reorder_levels(group_by_terms)

    # get the multi-index as a 2D array
    multi_index_array = np.vstack(multi_index.values)
    # create a string array of cell types for rollup
    # cell types are suffixed with ";;" to avoid collisions with cell type descendants
    # outside of each cell type's group.
    cell_types_for_rollup = _concatenate_columns_into_str(multi_index_array)

    # roll up the array along the first dimension (cell types;;groups)
    dense_tables = rollup_across_cell_type_descendants_array(dense_tables, cell_types_for_rollup)

    # sum, nnz, and sqsum come from the same sparse structure and are therefore guaranteed to be nonzero
    # in the same locations. We can use any of them to get the row and column coordinates of nonzero values.
    # rows and columns correspond to cell groups and genes respectively.
    row, col = dense_tables[:, :, 0].nonzero()

    new_df = pd.DataFrame(index=multi_index[row], data=genes[col], columns=["gene_ontology_term_id"]).reset_index()
    for i, col_name in enumerate(numeric_columns):
        new_df[col_name] = dense_tables[row, col, i]

    # get the total number of cells per tissue to populate the n_cells_tissue in the added entries
    n_cells_tissue = gene_expression_df.groupby("tissue_ontology_term_id").first()["n_cells_tissue"]

    new_df["n_cells_tissue"] = n_cells_tissue[new_df["tissue_ontology_term_id"]].values
    return new_df[gene_expression_df.columns]


def _concatenate_columns_into_str(value_arr):
    concatenated = np.char.add(value_arr[:, 0], ";;")
    for col in range(1, value_arr.shape[1]):
        if col > 1:
            concatenated = np.char.add(concatenated, "--")
        concatenated = np.char.add(concatenated, value_arr[:, col])
    return concatenated


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


@tracer.wrap(name="_rollup_cell_counts", service="wmg-api", resource="rollup", span_type="wmg-api")
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

    Redundant nodes are removed from the rolled up cell counts dataframe. A redundant node is a
    cell type that has no cell count values associated with it, but has ONLY ONE descendant cell type
    that has cell count values associated with it.

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
    if cell_counts_grouped_df.shape[0] == 0:
        return cell_counts_grouped_df

    rolled_up_cell_counts_grouped_df = cell_counts_grouped_df

    universal_set_cell_counts_grouped_df = _build_cell_count_groups_universal_set(cell_counts_grouped_df)
    index_names = universal_set_cell_counts_grouped_df.index.names

    universal_set_cell_counts_grouped_df = universal_set_cell_counts_grouped_df.reset_index()

    # rollup cell counts across cell type descendants
    rolled_up_cell_counts_grouped_df = rollup_across_cell_type_descendants(universal_set_cell_counts_grouped_df)

    rolled_up_cell_counts_grouped_df = rolled_up_cell_counts_grouped_df[
        rolled_up_cell_counts_grouped_df["n_cells_cell_type"] > 0
    ]

    rolled_up_cell_counts_grouped_df.set_index(index_names, inplace=True)

    index_names = list(index_names)
    index_names.remove("cell_type_ontology_term_id")
    index_names = ["cell_type_ontology_term_id"] + index_names
    multi_index = rolled_up_cell_counts_grouped_df.index.reorder_levels(index_names)

    cell_type_groups = _concatenate_columns_into_str(np.vstack(multi_index.values))

    cell_counts = dict(zip(cell_type_groups, rolled_up_cell_counts_grouped_df["n_cells_cell_type"].values))

    rolled_up_cell_counts_grouped_df = rolled_up_cell_counts_grouped_df[
        are_cell_types_not_redundant_nodes(cell_type_groups, cell_counts)
    ]
    return rolled_up_cell_counts_grouped_df
