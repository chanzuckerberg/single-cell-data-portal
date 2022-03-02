import logging
import os
import sys
from datetime import time
from typing import Dict

import tiledb

from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.wmg.cube_creation.loader import load
from backend.wmg.cube_creation.wmg_cube import create_cube

logger = logging.getLogger(__name__)

# TODO - does writing and reading directly from s3 slow down compute? test


def get_s3_uris():
    with db_session_manager() as session:
        dataset_ids = Dataset.list_ids_for_cube(session)
        s3_uris = DatasetAsset.list_s3_uris_for_datasets(session, dataset_ids, DatasetArtifactFileType.H5AD)
    return s3_uris


def copy_datasets_to_instance(dataset_directory):
    s3_uris = get_s3_uris()
    for dataset in s3_uris.keys():
        sync_command = f"aws s3 sync {s3_uris[dataset]} ./{dataset_directory}/{dataset}/local.h5ad"
        os.subprocess(sync_command) # TODO parallelize this step


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


def load_data_and_create_cube():
    copy_datasets_to_instance('wmg-datasets')
    load_datasets_into_corpus('wmg-datasets/', 'wmg-group')
    try:
        create_cube("wmg-group")
    except Exception as e:
        logger.error(f"Issue creating the cube: {e}")
    cell_type_by_tissue = get_cells_by_tissue_type("wmg-group")
    generate_cell_ordering(cell_type_by_tissue)
    update_s3_resources()


if __name__ == "__main__":
    rv = load_data_and_create_cube()
    sys.exit(rv)
