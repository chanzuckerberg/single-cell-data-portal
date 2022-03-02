import logging
import os
from datetime import time

from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.wmg.cube_creation.loader import load
logger = logging.getLogger(__name__)

# TODO - does writing and reading directly from s3 slow down compute? test


def get_s3_uris():
    with db_session_manager() as session:
        dataset_ids = Dataset.list_ids_for_cube(session)
        s3_uris = DatasetAsset.list_s3_uris_for_datasets(session, dataset_ids, DatasetArtifactFileType.H5AD)
    return s3_uris


def copy_datasets_to_instance():
    s3_uris = get_s3_uris()
    for dataset in s3_uris.keys():
        sync_command = f"aws s3 sync {s3_uris[dataset]} ./wmg-datasets/{dataset}/local.h5ad"
        os.subprocess(sync_command) # TODO parallelize this step


def load_datasets_into_corpus():
    try:
        load('wmg-datasets/', "wmg-group", True)
    except Exception as e:
        logger.error(f"Issue loading datasets into corpus: {e}")

def create_cube():
    pass

def generate_cell_ordering():
    pass


def update_cube():
    time_stamp = time.time()
    # create timestamp
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