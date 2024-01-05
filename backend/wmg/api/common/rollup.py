"""This module contains the implementation of the gene expression and cell count rollup feature.

The API public methods call the public methods in this module to perform the rollup operations.
"""


import numpy as np
import pandas as pd
from ddtrace import tracer
from pandas import DataFrame

from backend.common.utils.rollup import are_cell_types_not_redundant_nodes

######################### PUBLIC FUNCTIONS IN ALPHABETIC ORDER ##################################


@tracer.wrap(name="rollup", service="wmg-api", resource="query", span_type="wmg-api")
def rollup(df: pd.DataFrame, cell_type_ancestors: pd.Series, filter_redundant_nodes=True) -> DataFrame:
    if df.shape[0] == 0:
        return df

    df = df.copy()
    # cell type ontology term ids can be in the index or in a column
    cell_types = (
        df.index.get_level_values("cell_type_ontology_term_id")
        if "cell_type_ontology_term_id" not in df.columns
        else df["cell_type_ontology_term_id"].values
    )
    df["cell_type_ontology_term_id_ancestors"] = cell_type_ancestors[cell_types].values

    is_multi_index = isinstance(df.index, pd.MultiIndex)
    if is_multi_index:
        list(df.index.names)
        df = df.reset_index()

    del df["cell_type_ontology_term_id"]

    df = df.explode("cell_type_ontology_term_id_ancestors")

    df.rename(columns={"cell_type_ontology_term_id_ancestors": "cell_type_ontology_term_id"}, inplace=True)

    dim_cols = [col for col in df.columns if not np.issubdtype(df[col].dtype, np.number)]
    agg_dict = {col: "sum" if col != "n_cells_tissue" else "first" for col in df.columns if col not in dim_cols}
    rolled_up_df = df.groupby(dim_cols).agg(agg_dict)
    if filter_redundant_nodes:
        rolled_up_df = _filter_out_redundant_nodes(rolled_up_df)

    if not is_multi_index:
        rolled_up_df = rolled_up_df.reset_index()

    return rolled_up_df


######################### PRIVATE FUNCTIONS IN ALPHABETIC ORDER ##################################


def _filter_out_redundant_nodes(df):
    index_names = list(df.index.names)
    index_names.remove("cell_type_ontology_term_id")
    index_names = ["cell_type_ontology_term_id"] + index_names

    cell_type_keys = df.index.get_level_values(index_names[0]) + ";;"
    for i in index_names[1:]:
        cell_type_keys += df.index.get_level_values(i)

    cell_counts_dict = dict(zip(cell_type_keys, df["n_cells_cell_type"]))
    df = df[are_cell_types_not_redundant_nodes(cell_type_keys, cell_counts_dict)]
    return df
