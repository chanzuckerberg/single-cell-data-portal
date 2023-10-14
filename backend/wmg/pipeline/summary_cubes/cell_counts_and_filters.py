import json
import logging
import os

import numpy as np
import pandas as pd
import tiledb
from tiledbsoma import ExperimentAxisQuery

from backend.wmg.data.schemas.cube_schema import (
    cell_counts_logical_dims,
    cell_counts_schema,
)
from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    FILTER_RELATIONSHIPS_FILENAME,
)
from backend.wmg.data.utils import create_empty_cube, log_func_runtime, to_dict
from backend.wmg.pipeline.summary_cubes.constants import (
    CELL_COUNTS_CUBE_CREATED_FLAG,
    FILTER_RELATIONSHIPS_CREATED_FLAG,
)
from backend.wmg.pipeline.summary_cubes.utils import (
    load_pipeline_state,
    remove_accents,
    return_dataset_dict_w_publications,
    write_pipeline_state,
)

logger = logging.getLogger(__name__)


@log_func_runtime
def create_cell_counts_cube_and_filter_relationships(*, query: ExperimentAxisQuery, corpus_path: str):
    """
    Create cell count cube and write to disk
    """
    pipeline_state = load_pipeline_state(corpus_path)

    obs_df = query.obs().concat().to_pandas()

    logger.info("Creating the cell counts cube and filter relationships graph.")

    df = (
        obs_df.groupby(
            by=[dim for dim in cell_counts_logical_dims if dim != "publication_citation"],
            as_index=False,
        ).size()
    ).rename(columns={"size": "n_cells"})

    dataset_dict = return_dataset_dict_w_publications()
    df["publication_citation"] = [
        remove_accents(dataset_dict.get(dataset_id, "No Publication")) for dataset_id in df["dataset_id"]
    ]
    n_cells = df["n_cells"].to_numpy()
    df["n_cells"] = n_cells

    logger.info("Creating and writing filter relationships graph.")
    filter_relationships_linked_list = create_filter_relationships_graph(df)
    with open(f"{corpus_path}/{FILTER_RELATIONSHIPS_FILENAME}", "w") as f:
        json.dump(filter_relationships_linked_list, f)
    pipeline_state[FILTER_RELATIONSHIPS_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)

    uri = os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)
    create_empty_cube(uri, cell_counts_schema)
    logger.info("Writing cell counts cube.")
    tiledb.from_pandas(uri, df, mode="append")

    pipeline_state[CELL_COUNTS_CUBE_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)


@log_func_runtime
def create_filter_relationships_graph(df: pd.DataFrame) -> dict:
    """
    Create a graph of filter relationships

    Arguments
    ---------
    df: pd.DataFrame
        Dataframe containing the filter combinations per cell

    Returns
    -------
    filter_relationships_hash_table: dict
        A dictionary containing the filter relationships for each unique filter value.
        Filter values are their column name + __ + their value (e.g. "dataset_id__Single cell transcriptome analysis of human pancreas").
        The dictionary values are lists of filter values that are related to the key filter value. Relatedness indicates that these filters
        are co-occuring in at least one cell.
    """
    # get a dataframe of the columns that are not numeric
    df_filters = df.select_dtypes(exclude="number")
    # get a numpy array of the column names with shape (1, n_cols)
    cols = df_filters.columns.values[None, :]

    # tile the column names row-wise to match the shape of the dataframe and concatenate to the values
    # this ensures that filter values will never collide across columns.
    mat = np.tile(cols, (df.shape[0], 1)) + "__" + df_filters.values

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

    return filter_relationships_linked_list
