#!/usr/bin/env python3
import logging
import os
import sys

from backend.corpora.dataset_processing.process import (
    check_env,
    log_batch_environment,
    update_db,
    download_from_dropbox_url,
    validate_h5ad_file_and_add_labels,
    extract_metadata,
    cancel_dataset,
    notify_slack_failure,
    create_artifact,
    get_bucket_prefix,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

from backend.corpora.common.corpora_orm import (
    ProcessingStatus,
    DatasetArtifactFileType,
)
from backend.corpora.dataset_processing.exceptions import ProcessingCancelled, ProcessingFailed, ValidationFailed


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
    file_with_labels = validate_h5ad_file_and_add_labels(dataset_id, local_filename)
    # Process metadata
    metadata = extract_metadata(file_with_labels)
    update_db(dataset_id, metadata)

    # Upload the labeled dataset to the artifact bucket
    bucket_prefix = get_bucket_prefix(dataset_id)
    create_artifact(
        file_with_labels, DatasetArtifactFileType.H5AD, bucket_prefix, dataset_id, artifact_bucket, "h5ad_status"
    )


def main():
    return_value = 0
    check_env()
    log_batch_environment()

    dataset_id = os.environ["DATASET_ID"]
    dropbox_url = os.environ["DROPBOX_URL"]
    artifact_bucket = os.environ["ARTIFACT_BUCKET"]

    try:
        process(dataset_id, dropbox_url, artifact_bucket)
    except ProcessingCancelled:
        cancel_dataset(dataset_id)
    except (ValidationFailed, ProcessingFailed):
        logger.exception("An Error occurred while processing.")
        return_value = 1
    except Exception:
        message = f"Validation and labeling for dataset {dataset_id} failed."
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
