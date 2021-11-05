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
    process_cxg,
    get_bucket_prefix,
    LABELED_H5AD_FILENAME,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

from backend.corpora.common.corpora_orm import (
    ProcessingStatus,
)
from backend.corpora.dataset_processing.exceptions import ProcessingCancelled, ProcessingFailed, ValidationFailed


def download_from_s3(bucket_name: str, object_key: str, local_filename: str):
    s3_client.download_file(bucket_name, object_key, local_filename)


def process(dataset_id: str, artifact_bucket: str, labeled_h5ad_filename: str, cellxgene_bucket: str):
    """
    1. Download the labeled dataset from the artifact bucket
    1. Convert the labeled dataset to CXG
    2. Upload the CXG to the cellxgene bucket
    :param dataset_id:
    :param artifact_bucket:
    :param labeled_h5ad_filename:
    :param cellxgene_bucket:
    :return:
    """

    # Download the labeled dataset from the artifact bucket
    bucket_prefix = get_bucket_prefix(dataset_id)
    object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
    download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

    # Convert the labeled dataset to CXG and upload it to the cellxgene bucket
    process_cxg(labeled_h5ad_filename, dataset_id, cellxgene_bucket)


def main():
    return_value = 0
    check_env()
    log_batch_environment()

    dataset_id = os.environ["DATASET_ID"]
    artifact_bucket = os.environ["ARTIFACT_BUCKET"]
    cellxgene_bucket = os.environ["CELLXGENE_BUCKET"]

    try:
        process(dataset_id, artifact_bucket, LABELED_H5AD_FILENAME, cellxgene_bucket)
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
