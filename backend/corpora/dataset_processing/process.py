#!/usr/bin/env python3


"""
# Processing Status
## Initial
The initial processing_status when the container first runs is:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.WAITING
    upload_progress = 0
    upload_message = ""
    validation_status = ValidationStatus
    validation_message = ""
    rds_status = ConversionStatus
    cxg_status = ConversionStatus
    h5ad_status = ConversionStatus
}

## Upload

While uploading, upload_status changes UploadStatus.UPLOADING and upload_progress is updated regularly.
The processing_status should look like this:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADING
    upload_progress = 0.25
}

If upload succeeds the processing_status changes to:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
}


If upload fails the processing_status changes to:
{
    processing_status = ProcessingStatus.FAILURE
    upload_status = UploadStatus.FAILED
    upload_progress = 0.25
    upload_message = "Some message"
}

## Validation
After upload, validation starts and processing status changes to:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.VALIDATING
}

If validation succeeds the process_status changes to:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.VALID
    rds_status = ConversionStatus.CONVERTING
    cxg_status = ConversionStatus.CONVERTING
    h5ad_status = ConversionStatus.CONVERTING
}

If validation fails the processing_status change to:
{
    processing_status = ProcessingStatus.FAILURE
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.FAILED,
    h5ad_status = ConversionStatus.CONVERTED
}

## Conversion
As each conversion is started the status for the file format is set to CONVERTING then CONVERTED if the conversion
completes successfully. The status is then set to UPLOADING while the file is copied to s3 and UPLOADED on success.
 Cellxgene data is converted/uploaded first. The h5ad_status is set to CONVERTED after validation/label writing is
 complete.
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus
    rds_status = ConversionStatus.CONVERTING
    cxg_status = ConversionStatus.UPLOADED
    h5ad_status = ConversionStatus.CONVERTED
}

If a conversion fails the processing_status will indicated it as follow:
{
    processing_status = ProcessingStatus.PENDING
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.VALID
    rds_status = ConversionStatus.FAILED
    cxg_status = ConversionStatus.UPLOADED
    h5ad_status = ConversionStatus.UPLOADED
}

Once all conversion are complete, the conversion status for each file will be either UPLOADED or FAILED:
{
    processing_status = ProcessingStatus.FAILURE
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.VALID
    rds_status = ConversionStatus.FAILED
    cxg_status = ConversionStatus.UPLOADED
    h5ad_status = ConversionStatus.UPLOADED
}

# If upload, validation, and all conversions succeed, the overall dataset processing status will be set to SUCCESS:
{
    processing_status = ProcessingStatus.SUCCESS
    upload_status = UploadStatus.UPLOADED
    upload_progress = 1.0
    validation_status = ValidationStatus.VALID
    rds_status = ConversionStatus.UPLOADED
    cxg_status = ConversionStatus.UPLOADED
    h5ad_status = ConversionStatus.UPLOADED
}

# Standalone processing steps

## Seurat
This is used to recompute the Seurat artifact in place, starting from the original h5ad.
This is a state machine with a single state that mimics the Conversion step
of the main step function.

## CXG_Remaster
This is used to migrate the cxg to a different, more performant format. This is a state machine with a single
state that mimics the Conversion step of the main step function.

"""

import os
import sys

from backend.corpora.common.corpora_orm import (
    ConversionStatus,
    UploadStatus,
)
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.dataset_processing.common import update_db
from backend.corpora.dataset_processing.exceptions import (
    ProcessingCancelled,
    ProcessingFailed,
    ValidationFailed,
    ConversionFailed,
)
from backend.corpora.dataset_processing.logger import logger


def cancel_dataset(dataset_id):
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)
        dataset.asset_deletion()
        dataset.delete()
        logger.info("Upload Canceled.")


def log_batch_environment():
    batch_environment_variables = [
        "AWS_BATCH_CE_NAME",
        "AWS_BATCH_JOB_ATTEMPT",
        "AWS_BATCH_JOB_ID",
        "STEP_NAME",
        "DROPBOX_URL",  # TODO: Change to SOURCE_URI
        "ARTIFACT_BUCKET",
        "CELLXGENE_BUCKET",
        "DATASET_ID",
        "DEPLOYMENT_STAGE",
        "MAX_ATTEMPTS",
    ]
    env_vars = dict()
    for var in batch_environment_variables:
        env_vars[var] = os.getenv(var)
    logger.info(f"Batch Job Info: {env_vars}")


def main():
    log_batch_environment()
    dataset_id = os.environ["DATASET_ID"]
    step_name = os.environ["STEP_NAME"]
    return_value = 0
    logger.info(f"Processing dataset {dataset_id}")
    try:
        if step_name == "download-validate":
            from backend.corpora.dataset_processing.process_download_validate import (
                process,
            )

            process(dataset_id, os.environ["DROPBOX_URL"], os.environ["ARTIFACT_BUCKET"])
        elif step_name == "cxg":
            from backend.corpora.dataset_processing.process_cxg import process

            process(
                dataset_id,
                os.environ["ARTIFACT_BUCKET"],
                os.environ["CELLXGENE_BUCKET"],
            )
        elif step_name == "seurat":
            from backend.corpora.dataset_processing.process_seurat import process

            process(dataset_id, os.environ["ARTIFACT_BUCKET"])
        elif step_name == "cxg_remaster":
            from backend.corpora.dataset_processing.remaster_cxg import process

            process(dataset_id, os.environ["CELLXGENE_BUCKET"], dry_run=False)
        else:
            logger.error(f"Step function configuration error: Unexpected STEP_NAME '{step_name}'")

    except ProcessingCancelled:
        cancel_dataset(dataset_id)
    except (ValidationFailed, ProcessingFailed, ConversionFailed) as e:
        (status,) = e.args
        update_db(dataset_id, processing_status=status)
        logger.exception("An Error occurred while processing.")
        return_value = 1
    except Exception as e:
        logger.exception(f"An unexpected error occurred while processing the data set: {e}")
        (status,) = e.args
        if isinstance(status, dict):
            update_db(dataset_id, processing_status=status)
        else:
            if step_name == "download-validate":
                update_db(
                    dataset_id,
                    processing_status={"upload_status": UploadStatus.FAILED, "upload_message": str(e)},
                )
            elif step_name == "seurat":
                update_db(dataset_id, processing_status={"rds_status": ConversionStatus.FAILED})
            elif step_name == "cxg":
                update_db(dataset_id, processing_status={"cxg_status": ConversionStatus.FAILED})
        return_value = 1

    return return_value


if __name__ == "__main__":
    rv = main()
    sys.exit(rv)
