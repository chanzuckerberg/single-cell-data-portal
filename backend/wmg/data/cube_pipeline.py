import gc
import logging
import os
import pathlib
import sys
import time
from typing import List
import tiledb

from backend.wmg.data import extract
from backend.wmg.data.load_cube import update_s3_resources
from backend.wmg.data.load_corpus import load_h5ad
from backend.wmg.data.schemas.corpus_schema import create_tdb
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.transform import get_cells_by_tissue_type, generate_cell_ordering
from backend.wmg.data.wmg_cube import create_cube

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def load(dataset_directory: List, group_name: str, validate: bool = False):
    """
    Given the path to a directory containing one or more h5ad files and a group name, call the h5ad loading function
    on all files, loading/concatenating the datasets together under the group name
    """
    with tiledb.scope_ctx(create_ctx()):
        for dataset in os.listdir(dataset_directory):
            file_path = f"{dataset_directory}/{dataset}/local.h5ad"
            load_h5ad(
                file_path, group_name, validate
            )  # TODO Can this be parallelized? need to be careful handling global indexes but tiledb has a lock I think
            gc.collect()

        logger.info("all loaded, now consolidating.")
        for arr_name in [f"{group_name}/{name}" for name in ["obs", "var", "raw", "X"]]:
            tiledb.consolidate(arr_name)
            tiledb.vacuum(arr_name)


def load_data_and_create_cube(path_to_datasets: str, corpus_name: str = "corpus_group", snapshot_path=None):
    """
    Function to copy H5AD datasets (from a preconfiugred s3 bucket) to the path given then,
    open, transform, normalize and concatenate them together as a tiledb object with a global gene index
    under the given corpus name.
    A cube of expression summary statistics (known as the expression summary cube) across genes (but queryable by
    the given dimensions/data) is then generated and stored under big-cube.
    ## TODO add function to get cell count totals
    A per-tissue mapping of cell ontologies is generated and the files are copied to s3 under a shared timestamp,.
    On success the least recent set of files are removed from s3
    """
    if not snapshot_path:
        timestamp = int(time.time())
        snapshot_path = f"{pathlib.Path().resolve()}/{timestamp}"
    corpus_path = f"{snapshot_path}/{corpus_name}"
    if not tiledb.VFS().is_dir(corpus_path):
        create_tdb(snapshot_path, corpus_name)
    s3_uris = extract.get_dataset_s3_uris()
    extract.copy_datasets_to_instance(s3_uris, path_to_datasets)
    logger.info("Copied datasets to instance")

    load(path_to_datasets, corpus_path, True)
    logger.info("Loaded datasets into corpus")
    create_cube(corpus_path)
    logger.info("Built expression summary cube")

    cell_type_by_tissue = get_cells_by_tissue_type(corpus_path)
    generate_cell_ordering(cell_type_by_tissue)
    logger.info("Generated cell ordering json file")
    update_s3_resources()
    logger.info("Copied snapshot to s3")


if __name__ == "__main__":
    load_data_and_create_cube()
    sys.exit()
