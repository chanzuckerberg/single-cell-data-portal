import logging
import time
from typing import List

import tiledb

from backend.wmg.data.schemas.corpus_schema import OBS_ARRAY_NAME


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


def get_all_dataset_ids(tdb_group: str) -> List[str]:
    with tiledb.open(f"{tdb_group}/{OBS_ARRAY_NAME}", "r") as obs:
        all_dataset_ids = obs.query(attrs=[], dims=["dataset_id"]).df[:].dataset_id.unique()
    all_dataset_ids.sort()
    return all_dataset_ids


@log_func_runtime
def get_expression_summary_cube_gene_count(tbd_group: str) -> int:
    with tiledb.open(tbd_group) as obs:
        gene_count = len(obs.query(dims=["gene_ontology_term_id"]).df[:].gene_ontology_term_id.unique())
    return gene_count


@log_func_runtime
def get_cell_count_cube_count(tbd_group: str) -> int:
    with tiledb.open(tbd_group) as obs:
        cell_count = obs.query(attrs=["n_cells"]).df[:].n_cells_rollup.sum()
    return cell_count


def create_empty_cube(uri: str, schema):
    """
    Create an empty cube with expected schema (dimensions and attributes) at given uri
    """
    tiledb.Array.create(uri, schema, overwrite=True)
