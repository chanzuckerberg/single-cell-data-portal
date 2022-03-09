import logging
import os
import sys
import time
from typing import Dict
import subprocess
from pathlib import Path

import tiledb

from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.wmg.cube_creation.corpus_schema import create_tdb
from backend.wmg.cube_creation.loader import load
from backend.wmg.cube_creation.wmg_cube import create_cube


# TODO - does writing and reading directly from s3 slow down compute? test

Log_Format = "%(levelname)s %(asctime)s - %(message)s"

logging.basicConfig(
                    stream = sys.stdout,
                    filemode = "w",
                    format = Log_Format,
                    level = logging.DEBUG)

logger = logging.getLogger(__name__)


def get_s3_uris():
    with db_session_manager() as session:
        # TODO: remove when testing is completed
        # dataset_ids = []
        dataset_ids = Dataset.list_ids_for_cube(session)
        s3_uris = DatasetAsset.s3_uris_for_datasets(session, dataset_ids, DatasetArtifactFileType.H5AD)
    return s3_uris


def mock_s3_uris():
    return {
        # "cdd5d5f3-7ee5-44ed-95ad-54cb92aa5203": "s3://env-rdev-artifacts/ebezzi-test/cdd5d5f3-7ee5-44ed-95ad-54cb92aa5203/local.h5ad"
        "605073c5-5eca-45d7-8be7-571bf44bf746": "s3://env-rdev-wmg/ebezzi-wmg/605073c5-5eca-45d7-8be7-571bf44bf746/local.h5ad"
        }


def copy_datasets_to_instance(dataset_directory):
    s3_uris = mock_s3_uris()
    for dataset in s3_uris.keys():
        # Path(f"./{dataset_directory}/{dataset}/").mkdir(parents=True, exist_ok=True)
        sync_command = ["aws", "s3", "cp", s3_uris[dataset], f"./{dataset_directory}/{dataset}/local.h5ad"]
        # sync_command = f"aws s3 sync {s3_uris[dataset]} ./{dataset_directory}/{dataset}/local.h5ad"
        subprocess.run(sync_command) # TODO parallelize this step


def load_datasets_into_corpus(path_to_datasets, group_name):
    try:
        load(path_to_datasets, group_name, True)
    except Exception as e:
        logger.error(f"Issue loading datasets into corpus: {e}")


def get_cells_by_tissue_type(tdb_group: str) -> Dict:
    with tiledb.open(f"{tdb_group}/obs", "r") as obs:
     cell_tissue_types = (
         obs.query(
             attrs=[], dims=["tissue_ontology_term_id", "cell_type_ontology_term_id"]
         )
         .df[:]
         .drop_duplicates()
         .sort_values(by="tissue_ontology_term_id")
     )
    unique_tissue_ontology_term_id = cell_tissue_types.tissue_ontology_term_id.unique()
    cell_type_by_tissue = {}
    for x in unique_tissue_ontology_term_id:
        cell_type_by_tissue[x] = cell_tissue_types.loc[cell_tissue_types["tissue_ontology_term_id"] == x, "cell_type_ontology_term_id"]

    return cell_type_by_tissue

def generate_cell_ordering(cell_type_by_tissue):
    ## TODO port code from Emanuele notebook
    pass


def update_s3_resources():
    time_stamp = time.time()
    # copy cell ordering
    # copy corpus
    # copy cube
    update_latest_snapshot(time_stamp)
    remove_oldest_datasets()
    pass


def remove_oldest_datasets():
    pass


def update_latest_snapshot(time_stamp):
    pass


def load_data_and_create_cube(path_to_datasets, group_name):
    if not tiledb.VFS().is_dir(group_name):
        create_tdb(group_name)
    copy_datasets_to_instance(path_to_datasets)
    load_datasets_into_corpus(path_to_datasets, group_name)
    try:
        create_cube(group_name)
    except Exception as e:
        logger.error(f"Issue creating the cube: {e}")
    # cell_type_by_tissue = get_cells_by_tissue_type(group_name)
    # generate_cell_ordering(cell_type_by_tissue)
    # update_s3_resources()
    print("Cube creation script - called")
    return True


if __name__ == "__main__":
    rv = load_data_and_create_cube("datasets", "ebezzi-test-cube-wmg")
    sys.exit(rv)
