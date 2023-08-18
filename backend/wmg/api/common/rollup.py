"""This module contains the implementation of the gene expression and cell count rollup feature.

The API public methods call the public methods in this module to perform the rollup operations.
"""

import itertools
from typing import Tuple

import numpy as np
import pandas as pd
from pandas import DataFrame

from backend.common.utils.rollup import rollup_across_cell_type_descendants

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


def _add_missing_combinations_to_gene_expression_df_for_rollup(
    gene_expression_df, universal_set_cell_counts_df
) -> DataFrame:
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
    # extract group-by terms and queried genes from the input dataframes
    # if a queried gene is not present in the input dot plot dataframe, we can safely
    # ignore it as it need not be rolled up anyway.
    group_by_terms = list(universal_set_cell_counts_df.index.names)
    genes = list(set(gene_expression_df["gene_ontology_term_id"]))

    # get the names of the numeric columns
    numeric_columns = list(
        gene_expression_df.columns[[np.issubdtype(dtype, np.number) for dtype in gene_expression_df.dtypes]]
    )

    # exclude n_cells_tissue as we do not wish to roll it up
    if "n_cells_tissue" in numeric_columns:
        numeric_columns.remove("n_cells_tissue")

    # get the total number of cells per tissue to populate the n_cells_tissue in the added entries
    n_cells_tissue_dict = gene_expression_df.groupby("tissue_ontology_term_id").first()["n_cells_tissue"].to_dict()

    # get the set of available combinations of group-by terms from the aggregated cell counts
    available_combinations = set(universal_set_cell_counts_df.index.values)

    # for each gene, get the set of available combinations of group-by terms from the input expression dataframe
    entries_to_add = []
    for gene in genes:
        gene_expression_df_per_gene = gene_expression_df[gene_expression_df["gene_ontology_term_id"] == gene]
        available_combinations_per_gene = set(zip(*gene_expression_df_per_gene[group_by_terms].values.T))

        # get the combinations that are missing in the input expression dataframe
        # these combinations have no data but can be rescued by the roll-up operation
        missing_combinations = available_combinations.difference(available_combinations_per_gene)
        for combo in missing_combinations:
            entry = {dim: combo[i] for i, dim in enumerate(group_by_terms)}

            # If a tissue, T1, DOES NOT have ANY of the QUERIED GENES expressed, then entirely
            # omit all combos that contain tissue T1 from the candidate list of combos
            # that should be considered for rollup. Otherwise, include combos containing T1
            # in the candidate list of combos for rollup.
            if entry["tissue_ontology_term_id"] in n_cells_tissue_dict:
                entry.update({col: 0 for col in numeric_columns})
                entry["n_cells_tissue"] = n_cells_tissue_dict[entry["tissue_ontology_term_id"]]
                entry["gene_ontology_term_id"] = gene
                entries_to_add.append(entry)

    # add the missing entries to the input expression dataframe
    gene_expression_with_missing_combos_df = pd.concat((gene_expression_df, pd.DataFrame(entries_to_add)), axis=0)

    return gene_expression_with_missing_combos_df


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

        # Add index columns to the universal_set_cell_counts_grouped_df so that these columns can be
        # DIRECTLY accessed in the dataframe during the rollup operation while traversing the cell type descendants
        for col in universal_set_cell_counts_grouped_df.index.names:
            universal_set_cell_counts_grouped_df[col] = universal_set_cell_counts_grouped_df.index.get_level_values(col)

        # rollup cell counts across cell type descendants
        rolled_up_cell_counts_grouped_df = rollup_across_cell_type_descendants(universal_set_cell_counts_grouped_df)

        rolled_up_cell_counts_grouped_df = rolled_up_cell_counts_grouped_df[
            rolled_up_cell_counts_grouped_df["n_cells_cell_type"] > 0
        ]
        # Remove columns that were added to the cell counts dataframe for the purpose of rollup.
        # This is make it congruent with the structure of the input cell counts dataframe
        rolled_up_cell_counts_grouped_df.drop(columns=rolled_up_cell_counts_grouped_df.index.names, inplace=True)

    return rolled_up_cell_counts_grouped_df


def _rollup_gene_expression(gene_expression_df, universal_set_cell_counts_df) -> DataFrame:
    """
    Roll up numeric values across cell type descendants in the input gene expression dataframe.

    Accumulate gene expression values up the cell type ontology ancestor paths for every
    (gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>)
    combination in the input gene expression dataframe. The compare dimension is an optional
    field and in the case that it is missing, then accumulation happens for every
    (gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id) combination.

    For example:

    If (G1, T1, C1) has gene expression values in the input expression dataframe,
    and the tissue labeled T1 has another cell labeled C2 where C2 is an ANCESTOR of C1,
    then the rolled up gene expression dataframe must include CUMULATIVE gene expression values for the
    combination (G1, T1, C2) by adding in the gene expression for (G1, T1, C1).

    Parameters
    ----------
    gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe containing the dimensions across which the numeric columns will be
        aggregated.

    universal_set_cell_counts_df : pandas DataFrame
        Multi-indexed cell counts dataframe that contains "The Universal Set Of GroupBy Ontology Terms".

    Returns
    -------
    rolled_up_gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe with the same columns as the input gene expression dataframe,
        and likely greater size than the input gene expression dataframe, but with the numeric
        columns aggregated across the cell type's descendants.
    """
    rolled_up_gene_expression_df = gene_expression_df

    if gene_expression_df.shape[0] > 0:
        # For each gene in the query, add missing combinations (tissue, cell type, compare dimension)
        # to the expression dataframe
        gene_expression_with_missing_combos_df = _add_missing_combinations_to_gene_expression_df_for_rollup(
            gene_expression_df, universal_set_cell_counts_df
        )

        # Roll up expression dataframe
        rolled_up_gene_expression_df = rollup_across_cell_type_descendants(
            gene_expression_with_missing_combos_df, ignore_cols=["n_cells_tissue"]
        )

        # Filter out the entries that were added to the dataframe that remain zero after roll-up
        rolled_up_gene_expression_df = rolled_up_gene_expression_df[rolled_up_gene_expression_df["sum"] > 0]

    return rolled_up_gene_expression_df
