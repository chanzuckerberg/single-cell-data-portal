"""This module contains the functions that build the gene expression and cell counts data structures.

The API public methods call the public methods in this module to construct the gene expression
and cell count data structures process and return to the client.
"""

from typing import List, Tuple
from ddtrace import tracer
from pandas import DataFrame

DEFAULT_GROUP_BY_TERMS = ["tissue_ontology_term_id", "cell_type_ontology_term_id"]

######################### PUBLIC FUNCTIONS IN ALPHABETICAL ORDER ##################################


def agg_cell_type_counts(cell_counts: DataFrame, group_by_terms: List[str] = None) -> DataFrame:
    # Aggregate cube data by tissue, cell type
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS
    cell_counts_cell_type_agg = cell_counts.groupby(group_by_terms, as_index=True).sum(numeric_only=True)
    cell_counts_cell_type_agg.rename(columns={"n_total_cells": "n_cells_cell_type"}, inplace=True)
    return cell_counts_cell_type_agg


def agg_tissue_counts(cell_counts: DataFrame) -> DataFrame:
    # Aggregate cube data by tissue
    cell_counts_tissue_agg = cell_counts.groupby(["tissue_ontology_term_id"], as_index=True).sum(numeric_only=True)
    cell_counts_tissue_agg.rename(columns={"n_total_cells": "n_cells_tissue"}, inplace=True)
    return cell_counts_tissue_agg


def build_dot_plot_matrix(
    raw_gene_expression: DataFrame,
    cell_counts_cell_type_agg: DataFrame,
    cell_counts_tissue_agg: DataFrame,
    group_by_terms: List[str] = None,
) -> DataFrame:
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS

    # Aggregate cube data by gene, tissue, cell type
    expr_summary_agg = raw_gene_expression.groupby(["gene_ontology_term_id"] + group_by_terms, as_index=False).sum(
        numeric_only=True
    )
    return expr_summary_agg.join(cell_counts_cell_type_agg, on=group_by_terms, how="left").join(
        cell_counts_tissue_agg, on=["tissue_ontology_term_id"], how="left"
    )

@tracer.wrap(name="get_dot_plot_data", service="wmg-api", resource="query", span_type="wmg-api")
def get_dot_plot_data(
    raw_gene_expression: DataFrame,
    cell_counts: DataFrame,
    group_by_terms: List[str] = None,
) -> Tuple[DataFrame, DataFrame]:
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS
    # Get the dot plot matrix dataframe and aggregated cell counts per cell type
    cell_counts_cell_type_agg = agg_cell_type_counts(cell_counts, group_by_terms)
    cell_counts_tissue_agg = agg_tissue_counts(cell_counts)
    dot_plot_matrix_df = build_dot_plot_matrix(
        raw_gene_expression, cell_counts_cell_type_agg, cell_counts_tissue_agg, group_by_terms
    )
    return dot_plot_matrix_df, cell_counts_cell_type_agg
