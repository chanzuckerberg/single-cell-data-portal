"""This module contains the implementation of the gene expression and cell count rollup feature.

The API public methods call the public methods in this module to perform the rollup operations.
"""


import numpy as np
from ddtrace import tracer
from pandas import DataFrame

######################### PUBLIC FUNCTIONS IN ALPHABETIC ORDER ##################################


@tracer.wrap(name="rollup", service="wmg-api", resource="query", span_type="wmg-api")
def rollup(gene_expression_df, universal_set_cell_counts_df) -> DataFrame:
    # if the gene expression dataframe is empty, then there is nothing to roll up
    if gene_expression_df.shape[0] == 0 or universal_set_cell_counts_df.shape[0] == 0:
        return gene_expression_df, universal_set_cell_counts_df

    gene_expression_df["cell_type_ontology_term_id_ancestors"] = gene_expression_df.apply(lambda x: x.split(","))
    universal_set_cell_counts_df["cell_type_ontology_term_id_ancestors"] = universal_set_cell_counts_df.apply(
        lambda x: x.split(",")
    )

    del gene_expression_df["cell_type_ontology_term_id"]
    del universal_set_cell_counts_df["cell_type_ontology_term_id"]

    gene_expression_df = gene_expression_df.explode("cell_type_ontology_term_id_ancestors")
    universal_set_cell_counts_df = universal_set_cell_counts_df.explode("cell_type_ontology_term_id_ancestors")

    gene_expression_df.rename(
        columns={"cell_type_ontology_term_id_ancestors": "cell_type_ontology_term_id"}, inplace=True
    )
    universal_set_cell_counts_df.rename(
        columns={"cell_type_ontology_term_id_ancestors": "cell_type_ontology_term_id"}, inplace=True
    )
    # get non-numeric column names
    ge_dim_cols = [
        col for col in gene_expression_df.columns if not np.issubtype(gene_expression_df[col].dtype, np.number)
    ]
    cc_dim_cols = [
        col
        for col in universal_set_cell_counts_df.columns
        if not np.issubtype(universal_set_cell_counts_df[col].dtype, np.number)
    ]

    rolled_up_gene_expression_df = gene_expression_df.groupby(ge_dim_cols, as_index=False).sum()
    rolled_up_universal_set_cell_counts_df = universal_set_cell_counts_df.groupby(cc_dim_cols, as_index=False).sum()

    return rolled_up_gene_expression_df, rolled_up_universal_set_cell_counts_df
