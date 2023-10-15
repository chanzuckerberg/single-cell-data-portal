import json
import logging
import os

import numpy as np
import tiledb

from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    FILTER_RELATIONSHIPS_FILENAME,
)
from backend.wmg.data.utils import log_func_runtime, to_dict
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    FILTER_RELATIONSHIPS_CREATED_FLAG,
)
from backend.wmg.pipeline.utils import (
    load_pipeline_state,
    write_pipeline_state,
)

logger = logging.getLogger(__name__)


@log_func_runtime
def create_filter_relationships_graph(*, corpus_path: str) -> dict:
    """
    Create a graph of filter relationships

    Arguments
    ---------
    cell_counts_df: pd.DataFrame
        Dataframe containing the filter combinations per cell

    Returns
    -------
    filter_relationships_hash_table: dict
        A dictionary containing the filter relationships for each unique filter value.
        Filter values are their column name + __ + their value (e.g. "dataset_id__Single cell transcriptome analysis of human pancreas").
        The dictionary values are lists of filter values that are related to the key filter value. Relatedness indicates that these filters
        are co-occuring in at least one cell.
    """
    pipeline_state = load_pipeline_state(corpus_path)
    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise ValueError("'cell_counts' array does not exist. Please run 'create_cell_counts_cube' first.")

    logger.info("Creating filter relationships graph.")
    with tiledb.open(os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)) as cc_cube:
        cell_counts_df = cc_cube.df[:]

    # get a dataframe of the columns that are not numeric
    df_filters = cell_counts_df.select_dtypes(exclude="number")
    # get a numpy array of the column names with shape (1, n_cols)
    cols = df_filters.columns.values[None, :]

    # tile the column names row-wise to match the shape of the dataframe and concatenate to the values
    # this ensures that filter values will never collide across columns.
    mat = np.tile(cols, (cell_counts_df.shape[0], 1)) + "__" + df_filters.values

    # for each cell, get all pairwise combinations of filters compresent in that cell
    # these are the edges of the filter relationships graph
    Xs = []
    Ys = []
    for i in range(mat.shape[0]):
        Xs.extend(np.repeat(mat[i], mat[i].size))
        Ys.extend(np.tile(mat[i], mat[i].size))

    # get all the unique combinations of filters
    Xs, Ys = np.unique(np.array((Xs, Ys)), axis=1)

    # exclude self-relationships
    filt = Xs != Ys
    Xs, Ys = Xs[filt], Ys[filt]

    # convert the edges to a linked list representation
    filter_relationships_linked_list = to_dict(Xs, Ys)

    # reorganize the linked list representation to a nested linked list representation
    # where the filter columns are separated into distinct dictionaries
    # e.g. instead of {"cell_type_ontology_term_id__beta cell": ["dataset_id__Single cell transcriptome analysis of human pancreas", "assay_ontology_term_id__assay_type", ...], ...},
    # it's now {"cell_type_ontology_term_id__beta cell": {"dataset_id": ["dataset_id__Single cell transcriptome analysis of human pancreas", ...], "assay_ontology_term_id": ["assay_ontology_term_id__assay_type", ...], ...}, ...}.
    # This structure is easier to parse by the `/query` endpoint.
    for k, v in filter_relationships_linked_list.items():
        filter_relationships_linked_list[k] = to_dict([x.split("__")[0] for x in v], v)

    with open(f"{corpus_path}/{FILTER_RELATIONSHIPS_FILENAME}", "w") as f:
        json.dump(filter_relationships_linked_list, f)
    pipeline_state[FILTER_RELATIONSHIPS_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
