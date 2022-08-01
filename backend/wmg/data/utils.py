import logging
import time
from functools import cache
from typing import List

import tiledb

from backend.wmg.data.schemas.corpus_schema import OBS_ARRAY_NAME


@cache
def get_all_dataset_ids(tdb_group: str) -> List[str]:
    with tiledb.open(f"{tdb_group}/{OBS_ARRAY_NAME}", "r") as obs:
        all_dataset_ids = obs.query(attrs=[], dims=["dataset_id"]).df[:].dataset_id.unique()
    all_dataset_ids.sort()
    return all_dataset_ids


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


def create_empty_cube(uri: str, schema):
    """
    Create an empty cube with expected schema (dimensions and attributes) at given uri
    """
    tiledb.Array.create(uri, schema, overwrite=True)
