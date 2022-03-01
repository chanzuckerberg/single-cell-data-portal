import os

from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_session import db_session_manager


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



def copy_corpus_to_s3():
    pass

def copy_cube_to_s3():
    pass

def update_latest_snapshot():
    pass