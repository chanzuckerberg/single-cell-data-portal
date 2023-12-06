"""This module contains the implementation of the gene expression and cell count rollup feature.

The API public methods call the public methods in this module to perform the rollup operations.
"""


import numpy as np
from ddtrace import tracer
from pandas import DataFrame

from backend.common.utils.rollup import are_cell_types_not_redundant_nodes

######################### PUBLIC FUNCTIONS IN ALPHABETIC ORDER ##################################


@tracer.wrap(name="rollup", service="wmg-api", resource="query", span_type="wmg-api")
def rollup(gene_expression_df, cell_counts_df) -> DataFrame:
    # if the gene expression dataframe is empty, then there is nothing to roll up
    if gene_expression_df.shape[0] == 0 or cell_counts_df.shape[0] == 0:
        return gene_expression_df, cell_counts_df

    gene_expression_df["cell_type_ontology_term_id_ancestors"] = gene_expression_df[
        "cell_type_ontology_term_id_ancestors"
    ].apply(lambda x: x.split(","))
    cell_counts_df["cell_type_ontology_term_id_ancestors"] = cell_counts_df[
        "cell_type_ontology_term_id_ancestors"
    ].apply(lambda x: x.split(","))
    cell_counts_df = cell_counts_df.reset_index()

    del gene_expression_df["cell_type_ontology_term_id"]
    del cell_counts_df["cell_type_ontology_term_id"]

    gene_expression_df = gene_expression_df.explode("cell_type_ontology_term_id_ancestors")
    cell_counts_df = cell_counts_df.explode("cell_type_ontology_term_id_ancestors")

    gene_expression_df.rename(
        columns={"cell_type_ontology_term_id_ancestors": "cell_type_ontology_term_id"}, inplace=True
    )
    cell_counts_df.rename(columns={"cell_type_ontology_term_id_ancestors": "cell_type_ontology_term_id"}, inplace=True)
    # get non-numeric column names
    ge_dim_cols = [
        col for col in gene_expression_df.columns if not np.issubdtype(gene_expression_df[col].dtype, np.number)
    ]
    cc_dim_cols = [col for col in cell_counts_df.columns if not np.issubdtype(cell_counts_df[col].dtype, np.number)]

    agg_dict = {
        col: "sum" if col != "n_cells_tissue" else "first"
        for col in gene_expression_df.columns
        if col not in ge_dim_cols
    }
    rolled_up_gene_expression_df = gene_expression_df.groupby(ge_dim_cols).agg(agg_dict)
    rolled_up_cell_counts_df = cell_counts_df.groupby(cc_dim_cols).sum()

    return filter_out_redundant_nodes(rolled_up_gene_expression_df), filter_out_redundant_nodes(
        rolled_up_cell_counts_df
    )


def filter_out_redundant_nodes(df):
    index_names = list(df.index.names)
    index_names.remove("cell_type_ontology_term_id")
    index_names = ["cell_type_ontology_term_id"] + index_names

    cell_type_keys = df.index.get_level_values(index_names[0]) + ";;"
    for i in index_names[1:]:
        cell_type_keys += df.index.get_level_values(i)

    cell_counts_dict = dict(zip(cell_type_keys, df["n_cells_cell_type"]))
    df = df[are_cell_types_not_redundant_nodes(cell_type_keys, cell_counts_dict)]
    df = df[df.index.get_level_values("cell_type_ontology_term_id") != "Thing"]
    return df
