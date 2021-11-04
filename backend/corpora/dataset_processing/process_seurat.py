#!/usr/bin/env python3
import logging
import os
import sys

from backend.corpora.common.utils.s3_buckets import s3_client
from process import (
    check_env,
    log_batch_environment,
    update_db,
    cancel_dataset,
    notify_slack_failure,
    convert_file_ignore_exceptions,
    make_seurat,
    create_artifact,
    get_bucket_prefix,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

from backend.corpora.common.corpora_orm import (
    DatasetArtifactFileType,
    ProcessingStatus,
)
from backend.corpora.dataset_processing.exceptions import ProcessingCancelled, ProcessingFailed, ValidationFailed


def download_from_s3(bucket_name: str, object_key: str, local_filename: str):
    s3_client.download_file(bucket_name, object_key, local_filename)


def process(artifact_bucket: str, dataset_id: str, labeled_h5ad_filename: str, local_filename: str):

    bucket_prefix = get_bucket_prefix(dataset_id)
    object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"

    download_from_s3(artifact_bucket, object_key, local_filename)

    seurat_filename = convert_file_ignore_exceptions(
        make_seurat, local_filename, "Failed to convert dataset to Seurat format.", dataset_id, "rds_status"
    )

    if seurat_filename:
        create_artifact(
            seurat_filename, DatasetArtifactFileType.RDS, bucket_prefix, dataset_id, artifact_bucket, "rds_status"
        )


def main():
    return_value = 0
    check_env()
    log_batch_environment()

    artifact_bucket = os.environ["ARTIFACT_BUCKET"]
    dataset_id = os.environ["DATASET_ID"]
    labeled_h5ad_filename = os.environ["LABELED_H5AD_FILENAME"]
    local_filename = "local.h5ad"

    try:
        process(artifact_bucket, dataset_id, labeled_h5ad_filename, local_filename)
    except ProcessingCancelled:
        cancel_dataset(dataset_id)
    except (ValidationFailed, ProcessingFailed):
        logger.exception("An Error occurred while processing.")
        return_value = 1
    except Exception:
        message = "An unexpected error occurred while processing the data set."
        logger.exception(message)
        update_db(
            dataset_id, processing_status=dict(processing_status=ProcessingStatus.FAILURE, upload_message=message)
        )
        return_value = 1

    if return_value > 0:
        notify_slack_failure(dataset_id)
    return return_value


if __name__ == "__main__":
    rv = main()
    sys.exit(rv)
