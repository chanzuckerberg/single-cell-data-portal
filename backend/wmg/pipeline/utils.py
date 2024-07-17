import json
import logging
import os
import time
import unicodedata

import tiledb
from tiledb import ArraySchema

from backend.wmg.pipeline.constants import (
    DATASET_METADATA_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG,
    FILTER_RELATIONSHIPS_CREATED_FLAG,
    MARKER_GENES_CUBE_CREATED_FLAG,
    PIPELINE_STATE_FILENAME,
    PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG,
)

logger = logging.getLogger(__name__)


def load_pipeline_state(corpus_path: str):
    state_file = os.path.join(corpus_path, PIPELINE_STATE_FILENAME)
    if os.path.exists(state_file):
        with open(state_file, "r") as f:
            state = json.load(f)
    else:
        state = {
            EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG: False,
            EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG: False,
            MARKER_GENES_CUBE_CREATED_FLAG: False,
            FILTER_RELATIONSHIPS_CREATED_FLAG: False,
            PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG: False,
            DATASET_METADATA_CREATED_FLAG: False,
        }
    return state


def write_pipeline_state(pipeline_state: dict, corpus_path: str):
    with open(os.path.join(corpus_path, PIPELINE_STATE_FILENAME), "w") as f:
        json.dump(pipeline_state, f)


def remove_accents(input_str):
    nfkd_form = unicodedata.normalize("NFKD", input_str)
    return "".join([c for c in nfkd_form if not unicodedata.combining(c)])


def create_empty_cube_if_needed(uri: str, schema: ArraySchema):
    if not os.path.exists(uri):
        logger.info(f"Creating empty cube at {uri} with schema: {schema}")
        tiledb.Array.create(uri, schema, overwrite=True)


def log_func_runtime(func):
    # This decorator function logs the execution time of the function object passed
    def wrap_func(*args, **kwargs):
        logger = logging.getLogger(func.__module__)
        start = time.perf_counter()
        result = func(*args, **kwargs)
        stop = time.perf_counter()
        logger.info(f"Function {func.__name__} executed in {(stop-start):.4f}s")
        return result

    return wrap_func
