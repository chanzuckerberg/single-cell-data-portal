#!/usr/bin/env python3
import logging

from backend.corpora.dataset_processing.process import (
    convert_file_ignore_exceptions,
    download_from_s3,
    make_seurat,
    create_artifact,
    get_bucket_prefix,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

from backend.corpora.common.corpora_orm import DatasetArtifactFileType


def process(dataset_id: str, artifact_bucket: str):
    """
    1. Download the labeled dataset from the artifact bucket
    2. Convert it to Seurat format
    3. Upload the Seurat file to the artifact bucket
    :param artifact_bucket:
    :param dataset_id:
    :param labeled_h5ad_filename:
    :param local_filename:
    :return:
    """

    labeled_h5ad_filename = "local.h5ad"

    bucket_prefix = get_bucket_prefix(dataset_id)
    object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
    download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

    seurat_filename = convert_file_ignore_exceptions(
        make_seurat, labeled_h5ad_filename, "Failed to convert dataset to Seurat format.", dataset_id, "rds_status"
    )

    if seurat_filename:
        create_artifact(
            seurat_filename, DatasetArtifactFileType.RDS, bucket_prefix, dataset_id, artifact_bucket, "rds_status"
        )
