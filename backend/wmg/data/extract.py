# TODO generalize
import os

from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_session import db_session_manager

# TODO make dict of ontology ids and human readable names, store as a const elsewhere (in rds for easy updates?)
included_assay_ontology_ids = ['EFO:0008722', 'EFO:0010010', 'EFO:0010550', 'EFO:0010961', 'EFO:0030002',
                               'EFO:0009901', 'EFO:0011025', 'EFO:0009899', 'EFO:0009900', 'EFO:0009922',
                               'EFO_0030003', 'EFO:0030004', 'EFO:0008919']


def get_s3_uris():
    """
    Retrieve list of s3 uris for datasets included in the wmg cube
    """
    with db_session_manager() as session:

        dataset_ids = []
        published_dataset_non_null_assays = session.query(
            Dataset.table.id,
            Dataset.table.assay
        ).filter(
            Dataset.table.assay != 'null',
            Dataset.table.published == 'TRUE',
            Dataset.table.is_primary_data == 'PRIMARY',
            Dataset.table.collection_visibility == 'PUBLIC',
            Dataset.table.tombstone == 'f'
        ).all()
        for dataset in published_dataset_non_null_assays:
            if dataset[1]['ontology_term_id'] in included_assay_ontology_ids:
                dataset_ids.append(dataset[0])

        s3_uris = DatasetAsset.s3_uris_for_datasets(session, dataset_ids, DatasetArtifactFileType.H5AD)
    return s3_uris.keys()


def copy_datasets_to_instance(s3_uris, dataset_directory):
    """Copy given list of s3 uris to the provided path"""
    for dataset in s3_uris:
        sync_command = f"aws s3 sync {s3_uris[dataset]} ./{dataset_directory}/{dataset}/local.h5ad"
        os.subprocess(sync_command)  # TODO parallelize this step
