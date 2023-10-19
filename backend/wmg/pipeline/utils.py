import json
import logging
import os
import time
import unicodedata

import numpy as np
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
from backend.wmg.pipeline.utils import (
    get_collections_from_discover_api,
    get_datasets_from_discover_api,
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


def return_dataset_dict_w_publications():
    datasets = get_datasets_from_discover_api()
    collections = get_collections_from_discover_api()
    collections_dict = {collection["collection_id"]: collection for collection in collections}

    # cchoi: creating a helper function to format citations properly
    def create_formatted_citation(collection):
        publisher_metadata = collection["publisher_metadata"]
        if publisher_metadata is None:
            return "No Publication"
        first_author = collection["publisher_metadata"]["authors"][0]
        # first_author could be either 'family' or 'name'
        citation = f"{first_author['family'] if 'family' in first_author else first_author['name']} et al. {collection['publisher_metadata']['journal']} {collection['publisher_metadata']['published_year']}"
        formatted_citation = "No Publication" if collection["publisher_metadata"]["is_preprint"] else citation
        return formatted_citation

    dataset_dict = {}
    for dataset in datasets:
        dataset_id = dataset["dataset_id"]
        collection = collections_dict[dataset["collection_id"]]
        dataset_dict[dataset_id] = create_formatted_citation(collection)

    return dataset_dict


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


def to_dict(a, b):
    """
    convert a flat key array (a) and a value array (b) into a dictionary with values grouped by keys
    """
    a = np.array(a)
    b = np.array(b)
    idx = np.argsort(a)
    a = a[idx]
    b = b[idx]
    bounds = np.where(a[:-1] != a[1:])[0] + 1
    bounds = np.append(np.append(0, bounds), a.size)
    bounds_left = bounds[:-1]
    bounds_right = bounds[1:]
    slists = [b[bounds_left[i] : bounds_right[i]] for i in range(bounds_left.size)]
    d = dict(zip(np.unique(a), [list(set(x)) for x in slists]))
    return d
