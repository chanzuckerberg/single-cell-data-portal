#!/usr/bin/env python3
import logging

from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.dataset_processing.process import (
    convert_file_ignore_exceptions,
    download_from_s3,
    make_seurat,
    create_artifact,
    get_bucket_prefix,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

from backend.corpora.common.corpora_orm import ConversionStatus, DatasetArtifactFileType


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

    # If the validator previously marked the dataset as rds_status.SKIPPED, do not start the Seurat processing
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)

        if dataset.processing_status.rds_status == ConversionStatus.SKIPPED:
            logger.info("Skipping Seurat conversion")
            return

    labeled_h5ad_filename = "local.h5ad"

    bucket_prefix = get_bucket_prefix(dataset_id)
    object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
    logger.warning(f"Download {object_key} from S3 bucket {artifact_bucket}")
    download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

    seurat_filename = convert_file_ignore_exceptions(
        make_seurat, labeled_h5ad_filename, "Failed to convert dataset to Seurat format.", dataset_id, "rds_status"
    )

    if seurat_filename:
        create_artifact(
            seurat_filename, DatasetArtifactFileType.RDS, bucket_prefix, dataset_id, artifact_bucket, "rds_status"
        )
