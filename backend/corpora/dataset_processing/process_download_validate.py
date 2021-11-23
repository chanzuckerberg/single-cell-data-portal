#!/usr/bin/env python3
import logging

from backend.corpora.dataset_processing.process import (
    update_db,
    download_from_dropbox_url,
    validate_h5ad_file_and_add_labels,
    extract_metadata,
    create_artifact,
    get_bucket_prefix,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

from backend.corpora.common.corpora_orm import (
    ProcessingStatus,
    DatasetArtifactFileType,
)


def process(dataset_id: str, dropbox_url: str, artifact_bucket: str):
    """
    1. Download the original dataset from Dropbox
    2. Validate and label it
    3. Upload the labeled dataset to the artifact bucket
    :param dataset_id:
    :param dropbox_url:
    :param artifact_bucket:
    :return:
    """

    update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.PENDING))

    # Download the original dataset from Dropbox
    local_filename = download_from_dropbox_url(
        dataset_id,
        dropbox_url,
        "raw.h5ad",
    )

    # TODO: We should store the original dataset file on S3 for provenance.
    #       This will allow for easy reconversion of the corpus in the near future.

    # Validate and label the dataset
    file_with_labels, can_convert_to_seurat = validate_h5ad_file_and_add_labels(dataset_id, local_filename)
    # Process metadata
    metadata = extract_metadata(file_with_labels)
    update_db(dataset_id, metadata)

    # Upload the labeled dataset to the artifact bucket
    bucket_prefix = get_bucket_prefix(dataset_id)
    create_artifact(
        file_with_labels, DatasetArtifactFileType.H5AD, bucket_prefix, dataset_id, artifact_bucket, "h5ad_status"
    )
