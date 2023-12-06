"""This module contains the implementation of the gene expression and cell count rollup feature.

The API public methods call the public methods in this module to perform the rollup operations.
"""


from ddtrace import tracer
from pandas import DataFrame

from backend.common.utils.rollup import (
    rollup_across_cell_type_descendants,
)

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

    rolled_up_gene_expression_df = rollup_across_cell_type_descendants(gene_expression_df)
    rolled_up_universal_set_cell_counts_df = rollup_across_cell_type_descendants(universal_set_cell_counts_df)

    return rolled_up_gene_expression_df, rolled_up_universal_set_cell_counts_df
