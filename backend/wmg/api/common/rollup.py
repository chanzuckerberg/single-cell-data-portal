"""This module contains the implementation of the gene expression and cell count rollup feature.

The API public methods call the public methods in this module to perform the rollup operations.
"""


from typing import Union

import numpy as np
import pandas as pd
from ddtrace import tracer
from pandas import DataFrame

from backend.common.utils.rollup import are_cell_types_not_redundant_nodes

######################### PUBLIC FUNCTIONS IN ALPHABETIC ORDER ##################################


@tracer.wrap(name="rollup", service="wmg-api", resource="query", span_type="wmg-api")
def rollup(df: pd.DataFrame, cell_type_ancestors: Union[pd.Series, None], filter_redundant_nodes=True) -> DataFrame:
    """
    This function performs rollup operations on a given DataFrame. It aggregates the data based on cell type ontology
    term ids and their ancestors. It also provides an option to filter out redundant nodes.
    Redundance means that a cell type has the same number of cells as one of its descendants. This indicates that
    the cell type contains exactly the same set of cells as one of its descendant subtrees. The only way this can occur is
    if the cell type belongs to a linear chain of cell types, all of which have the same number of cells. In these cases,
    the cell type is redundant and should be excluded. We only wish to keep the leaf nodes of these linear chains.

    Parameters
    -----------
    df : pd.DataFrame
        The input DataFrame on which rollup operation is to be performed.
    cell_type_ancestors : Union[pd.Series, None]
        A pandas Series containing cell type ancestors or None. None is expected for unit tests that are not testing for
        correctness of rollup and do not have the cell type ancestors artifact available.
    filter_redundant_nodes : bool, optional
        A flag to indicate whether to filter out redundant nodes, by default True.

    Returns
    --------
    DataFrame
        The rolled up DataFrame.
    """

    if df.shape[0] == 0:
        return df

    df = df.copy()
    # cell type ontology term ids can be in the index or in a column
    cell_types = (
        df.index.get_level_values("cell_type_ontology_term_id")
        if "cell_type_ontology_term_id" not in df.columns
        else df["cell_type_ontology_term_id"].values
    )
    df["cell_type_ontology_term_id_ancestors"] = (
        cell_type_ancestors[cell_types].values if cell_type_ancestors is not None else cell_types
    )

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
    """
    Filters out redundant nodes from the DataFrame.

    Parameters
    -----------
    df : DataFrame
        The DataFrame to filter.

    Returns
    --------
    DataFrame
        The filtered DataFrame.
    """

    index_names = list(df.index.names)
    index_names.remove("cell_type_ontology_term_id")
    index_names = ["cell_type_ontology_term_id"] + index_names

    cell_type_keys = df.index.get_level_values(index_names[0]) + ";;"
    for i in index_names[1:]:
        cell_type_keys += df.index.get_level_values(i)

    cell_counts_dict = dict(zip(cell_type_keys, df["n_cells_cell_type"]))
    df = df[are_cell_types_not_redundant_nodes(cell_type_keys, cell_counts_dict)]
    return df
