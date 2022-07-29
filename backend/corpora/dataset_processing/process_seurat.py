#!/usr/bin/env python3
import logging

from backend.corpora.common.entities.tiledb_data import TileDBData
from backend.corpora.dataset_processing.process import (
    convert_file_ignore_exceptions,
    download_from_s3,
    make_seurat,
    create_artifact,
    get_bucket_prefix,
    replace_artifact,
)
from datetime import datetime

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
    db = TileDBData()

    # If the validator previously marked the dataset as rds_status.SKIPPED, do not start the Seurat processing
    dataset = db.get_dataset(dataset_id)
    if dataset['processing_status']['rds_status'] == ConversionStatus.SKIPPED.name:
        logger.info("Skipping Seurat conversion")
        return

    labeled_h5ad_filename = "local.h5ad"

    # Four cases:
    # 1. newly processed dataset: h5ad will exist with key == dataset_id, rds will not exist
    # 2. non-revised reprocessed dataset with seurat
    # 3. revised reprocessed dataset with no seurat: h5ad will exist with key != dataset_id, rds will not exist
    # 4. revised reprocessed dataset with "to be replaced" seurat: h5ad will exist with key ≠ dataset_id,
    #    rds will exist

    h5ad_uri = next(a['s3_uri'] for a in dataset['dataset_assets'] if a['filetype'] == DatasetArtifactFileType.H5AD.name)
    rds_artifacts = [a for a in dataset['dataset_assets'] if a['filetype'] == DatasetArtifactFileType.RDS.name]
    artifact_key = h5ad_uri.split("/")[-2]

    if artifact_key == dataset_id and not rds_artifacts:  # Case 1

        bucket_prefix = get_bucket_prefix(dataset_id)
        object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
        download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

        seurat_filename = convert_file_ignore_exceptions(
            make_seurat,
            labeled_h5ad_filename,
            "Failed to convert dataset to Seurat format.",
            dataset_id,
            "rds_status",
        )

        if seurat_filename:
            create_artifact(
                seurat_filename,
                DatasetArtifactFileType.RDS.name,
                bucket_prefix,
                dataset_id,
                artifact_bucket,
                "rds_status",
            )

    elif artifact_key == dataset_id and rds_artifacts:  # Case 2
        logger.warning(f"Reprocessing Seurat for existing dataset {artifact_key}, will replace the S3 file only")
        bucket_prefix = get_bucket_prefix(dataset_id)
        object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
        download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

        seurat_filename = convert_file_ignore_exceptions(
            make_seurat,
            labeled_h5ad_filename,
            "Failed to convert dataset to Seurat format.",
            dataset_id,
            "rds_status",
        )

        if seurat_filename:
            replace_artifact(seurat_filename, bucket_prefix, artifact_bucket)
            rds_artifact = rds_artifacts[0]  # Only one RDS artifact for dataset will ever exist
            rds_artifact.updated_at = datetime.utcnow()

    elif artifact_key != dataset_id and not rds_artifacts:  # Case 3
        logger.warning(f"Found existing artifacts in {artifact_key} but no RDS, creating a new artifact")
        bucket_prefix = get_bucket_prefix(artifact_key)
        object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
        download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

        seurat_filename = convert_file_ignore_exceptions(
            make_seurat,
            labeled_h5ad_filename,
            "Failed to convert dataset to Seurat format.",
            dataset_id,
            "rds_status",
        )

        if seurat_filename:
            create_artifact(
                seurat_filename,
                DatasetArtifactFileType.RDS.name,
                bucket_prefix,
                dataset_id,
                artifact_bucket,
                "rds_status",
            )

    else:  # Case 4
        logger.warning(f"Found existing artifacts in {artifact_key} including an RDS, replacing it")

        bucket_prefix = get_bucket_prefix(artifact_key)
        object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
        download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

        seurat_filename = convert_file_ignore_exceptions(
            make_seurat,
            labeled_h5ad_filename,
            "Failed to convert dataset to Seurat format.",
            dataset_id,
            "rds_status",
        )

        if seurat_filename:
            replace_artifact(seurat_filename, bucket_prefix, artifact_bucket)
            rds_artifact = rds_artifacts[0]  # Only one RDS artifact for dataset will ever exist
            rds_artifact.updated_at = datetime.utcnow()
