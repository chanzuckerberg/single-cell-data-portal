import gc
import logging
import os
import pathlib
import sys
import time
from typing import List
import tiledb

from backend.wmg.data import extract
from backend.wmg.data.load_cube import update_s3_resources, upload_artifacts_to_s3
from backend.wmg.data.load_corpus import load_h5ad
from backend.wmg.data.schemas.corpus_schema import create_tdb, INTEGRATED_ARRAY_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.transform import (
    generate_primary_filter_dimensions,
    get_cell_types_by_tissue,
    generate_cell_ordering,
)
from backend.wmg.data.validation.validation import Validation
from backend.wmg.data.wmg_cube import create_cubes

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def load(dataset_directory: List, corpus_path: str, validate: bool = False):
    """
    Given the path to a directory containing one or more h5ad files and a group name, call the h5ad loading function
    on all files, loading/concatenating the datasets together under the group name
    """
    with tiledb.scope_ctx(create_ctx()):
        dataset_count = len(os.listdir(dataset_directory))
        i = 0
        for dataset in os.listdir(dataset_directory):
            i += 1
            logger.info(f"Processing dataset {i} of {dataset_count}")
            h5ad_file_path = f"{dataset_directory}/{dataset}/local.h5ad"
            load_h5ad(
                h5ad_file_path, corpus_path, validate
            )  # TODO Can this be parallelized? need to be careful handling global indexes but tiledb has a lock I think
            gc.collect()

        logger.info("all loaded, now consolidating.")
        for arr_name in [f"{corpus_path}/{name}" for name in ["obs", "var", INTEGRATED_ARRAY_NAME]]:
            tiledb.consolidate(arr_name)
            tiledb.vacuum(arr_name)


def load_data_and_create_cube(
    path_to_h5ad_datasets: str,
    corpus_name: str = "corpus_group",
    snapshot_path=None,
    extract_data=True,
    validate_cubes=True,
):
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

    if extract_data:
        s3_uris = extract.get_dataset_s3_uris()
        extract.copy_datasets_to_instance(s3_uris, path_to_h5ad_datasets)
        logger.info("Copied datasets to instance")

    load(path_to_h5ad_datasets, corpus_path, True)
    logger.info("Loaded datasets into corpus")
    create_cubes(corpus_path)
    logger.info("Built expression summary cube")
    if validate_cubes:
        is_valid = Validation(corpus_path).validate_cube()
        if is_valid is False:
            return
    cell_type_by_tissue = get_cell_types_by_tissue(corpus_path)
    generate_cell_ordering(snapshot_path, cell_type_by_tissue)
    generate_primary_filter_dimensions(snapshot_path, corpus_name, timestamp)
    logger.info("Generated cell ordering json file")
    upload_artifacts_to_s3(snapshot_path, timestamp)
    logger.info("Copied snapshot to s3")
    if validate_cubes:
        update_s3_resources(timestamp)
    logger.info(f"Updated latest_snapshot_identifier in s3. Current snapshot id is {timestamp}")


if __name__ == "__main__":
    load_data_and_create_cube("datasets", ".")
    sys.exit()
