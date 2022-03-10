import os

from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_session import db_session_manager


def get_dataset_s3_uris():
    """
    Retrieve list of s3 uris for datasets included in the wmg cube
    """
    with db_session_manager() as session:
        dataset_ids = Dataset.list_ids_for_cube(session)
        s3_uris = DatasetAsset.s3_uris_for_datasets(session, dataset_ids, DatasetArtifactFileType.H5AD)
    return s3_uris.keys()


def copy_datasets_to_instance(s3_uris, dataset_directory):
    """Copy given list of s3 uris to the provided path"""
    for dataset in s3_uris:
        sync_command = f"aws s3 sync {s3_uris[dataset]} ./{dataset_directory}/{dataset}/local.h5ad"
        os.subprocess(sync_command)  # TODO parallelize this step
