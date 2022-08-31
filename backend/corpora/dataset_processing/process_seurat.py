#!/usr/bin/env python3
import os
import subprocess

from backend.corpora.common.entities import DatasetAsset
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.dataset_processing.logger import logger
from backend.corpora.dataset_processing.common import create_artifact, download_from_s3, convert_file, get_bucket_prefix
from datetime import datetime


from backend.corpora.common.corpora_orm import ConversionStatus, DatasetArtifactFileType


def process(dataset_id: str, artifact_bucket: str):
    """
    1. Download the labeled dataset from the artifact bucket
    2. Convert it to Seurat format
    3. Upload the Seurat file to the artifact bucket
    :param artifact_bucket:
    :param dataset_id:
    :return:
    """

    # If the validator previously marked the dataset as rds_status.SKIPPED, do not start the Seurat processing
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)

        if dataset.processing_status.rds_status == ConversionStatus.SKIPPED:
            logger.info("Skipping Seurat conversion")
            return

        labeled_h5ad_filename = "local.h5ad"

        # Four cases:
        # 1. newly processed dataset: h5ad will exist with key == dataset_id, rds will not exist
        # 2. non-revised reprocessed dataset with seurat
        # 3. revised reprocessed dataset with no seurat: h5ad will exist with key != dataset_id, rds will not exist
        # 4. revised reprocessed dataset with "to be replaced" seurat: h5ad will exist with key ≠ dataset_id,
        #    rds will exist

        h5ad_uri = next(a.s3_uri for a in dataset.artifacts if a.filetype == DatasetArtifactFileType.H5AD)
        rds_artifacts = [a for a in dataset.artifacts if a.filetype == DatasetArtifactFileType.RDS]
        artifact_key = h5ad_uri.split("/")[-2]

        if artifact_key == dataset_id and not rds_artifacts:  # Case 1

            bucket_prefix = get_bucket_prefix(dataset_id)
            object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
            download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

            seurat_filename = convert_file(
                make_seurat,
                labeled_h5ad_filename,
                "Failed to convert dataset to Seurat format.",
                dataset_id,
                "rds_status",
            )

            if seurat_filename:
                create_artifact(
                    seurat_filename,
                    DatasetArtifactFileType.RDS,
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

            seurat_filename = convert_file(
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

            seurat_filename = convert_file(
                make_seurat,
                labeled_h5ad_filename,
                "Failed to convert dataset to Seurat format.",
                dataset_id,
                "rds_status",
            )

            if seurat_filename:
                create_artifact(
                    seurat_filename,
                    DatasetArtifactFileType.RDS,
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

            seurat_filename = convert_file(
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


def make_seurat(local_filename):
    """Create a Seurat rds file from the AnnData file."""

    try:
        subprocess.run(
            [
                "Rscript",
                os.path.join(os.path.abspath(os.path.dirname(__file__)), "make_seurat.R"),
                local_filename,
            ],
            capture_output=True,
            check=True,
        )
    except subprocess.CalledProcessError as ex:
        msg = f"Seurat conversion failed: {ex.output} {ex.stderr}"
        logger.exception(msg)
        raise RuntimeError(msg) from ex

    return local_filename.replace(".h5ad", ".rds")


def replace_artifact(
    file_name: str,
    bucket_prefix: str,
    artifact_bucket: str,
):
    logger.info(f"Uploading [{bucket_prefix}/{file_name}] to S3 bucket: [{artifact_bucket}].")
    DatasetAsset.upload(file_name, bucket_prefix, artifact_bucket)
