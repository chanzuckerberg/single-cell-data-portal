#!/usr/bin/env python3
import sys
import os
from process import (
    check_env,
    log_batch_environment,
    update_db,
    download_from_dropbox_url,
    validate_h5ad_file_and_add_labels,
    extract_metadata,
    cancel_dataset,
    notify_slack_failure,
    process_cxg,
    get_bucket_prefix,
    create_artifact
)
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.corpora_orm import (
    ConversionStatus,
    DatasetArtifactFileType,
    DatasetArtifactType,
    ProcessingStatus,
    ValidationStatus,
)
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_helpers import processing_status_updater
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.dl_sources.url import from_url
from backend.corpora.dataset_processing.download import download
from backend.corpora.dataset_processing.exceptions import ProcessingCancelled, ProcessingFailed, ValidationFailed
from backend.corpora.dataset_processing.h5ad_data_file import H5ADDataFile
from backend.corpora.dataset_processing.slack import format_slack_message


def process(dataset_id, dropbox_url, cellxgene_bucket, artifact_bucket):

    process_cxg(file_with_labels, dataset_id, cellxgene_bucket)
    bucket_prefix = get_bucket_prefix(dataset_id)
    logger.info("Creating Artifacts.")
    # upload AnnData
    create_artifact(
        local_filename, DatasetArtifactFileType.H5AD, bucket_prefix, dataset_id, artifact_bucket, "h5ad_status"
    )
    # update_db(dataset_id, processing_status=dict(processing_status=ProcessingStatus.SUCCESS))


def main():
    return_value = 0
    check_env()
    log_batch_environment()
    dataset_id = os.environ["DATASET_ID"]
    try:
        process(dataset_id, os.environ["DROPBOX_URL"], os.environ["CELLXGENE_BUCKET"], os.environ["ARTIFACT_BUCKET"])
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
