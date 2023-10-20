import json
import logging
import os

import tiledb

from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    FILTER_RELATIONSHIPS_FILENAME,
)
from backend.wmg.data.utils import build_filter_relationships
from backend.wmg.pipeline.constants import (
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    FILTER_RELATIONSHIPS_CREATED_FLAG,
)
from backend.wmg.pipeline.errors import PipelineStepMissing
from backend.wmg.pipeline.utils import load_pipeline_state, log_func_runtime, write_pipeline_state

logger = logging.getLogger(__name__)


@log_func_runtime
def create_filter_relationships_graph(corpus_path: str) -> dict:
    """
    Create a graph of filter relationships

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
        raise PipelineStepMissing("cell_counts")

    logger.info("Creating the filter relationships graph.")
    with tiledb.open(os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)) as cc_cube:
        cell_counts_df = cc_cube.df[:]

    filter_relationships_linked_list = build_filter_relationships(cell_counts_df)

    with open(f"{corpus_path}/{FILTER_RELATIONSHIPS_FILENAME}", "w") as f:
        json.dump(filter_relationships_linked_list, f)
    pipeline_state[FILTER_RELATIONSHIPS_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)
