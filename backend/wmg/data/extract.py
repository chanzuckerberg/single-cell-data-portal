import subprocess

from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import Collection, Dataset, DatasetAsset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.wmg.data.wmg_constants import included_assay_ontologies


def get_dataset_s3_uris():
    """
    Retrieve list of s3 uris for datasets included in the wmg cube
    """
    with db_session_manager() as session:

        dataset_ids = []
        published_dataset_non_null_assays = (
            session.query(Dataset.table.id, Dataset.table.assay)
            .join(Dataset.table.collection)
            .filter(
                Dataset.table.assay != "null",
                Dataset.table.published == "TRUE",
                Dataset.table.is_primary_data == "PRIMARY",
                Collection.table.visibility == "PUBLIC",
                Dataset.table.tombstone == "FALSE",
            )
            .all()
        )
        for dataset in published_dataset_non_null_assays:
            if dataset[1]["ontology_term_id"] in included_assay_ontologies:
                dataset_ids.append(dataset[0])

        s3_uris = DatasetAsset.s3_uris_for_datasets(session, dataset_ids, DatasetArtifactFileType.H5AD)
    return s3_uris


def copy_datasets_to_instance(s3_uris, dataset_directory):
    """Copy given list of s3 uris to the provided path"""
    for dataset in s3_uris:
        copy_command = ["aws", "s3", "cp", s3_uris[dataset], f"./{dataset_directory}/{dataset}/local.h5ad"]
        subprocess.run(copy_command)
