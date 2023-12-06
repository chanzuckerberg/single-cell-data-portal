"""This module contains the implementation of the gene expression and cell count rollup feature.

The API public methods call the public methods in this module to perform the rollup operations.
"""


import numpy as np
from ddtrace import tracer
from pandas import DataFrame

######################### PUBLIC FUNCTIONS IN ALPHABETIC ORDER ##################################


@tracer.wrap(name="rollup", service="wmg-api", resource="query", span_type="wmg-api")
def rollup(gene_expression_df, cell_counts_df) -> DataFrame:
    # if the gene expression dataframe is empty, then there is nothing to roll up
    if gene_expression_df.shape[0] == 0 or cell_counts_df.shape[0] == 0:
        return gene_expression_df, cell_counts_df

    gene_expression_df["cell_type_ontology_term_id_ancestors"] = gene_expression_df.apply(lambda x: x.split(","))
    cell_counts_df["cell_type_ontology_term_id_ancestors"] = cell_counts_df.apply(lambda x: x.split(","))

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
        col for col in gene_expression_df.columns if not np.issubtype(gene_expression_df[col].dtype, np.number)
    ]
    cc_dim_cols = [col for col in cell_counts_df.columns if not np.issubtype(cell_counts_df[col].dtype, np.number)]

    agg_dict = {
        col: "sum" if col != "n_cells_tissue" else "first"
        for col in gene_expression_df.columns
        if col not in ge_dim_cols
    }
    rolled_up_gene_expression_df = gene_expression_df.groupby(ge_dim_cols).agg(agg_dict)
    rolled_up_cell_counts_df = cell_counts_df.groupby(cc_dim_cols).sum()

    return rolled_up_gene_expression_df, rolled_up_cell_counts_df
