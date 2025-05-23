import json
import logging
import os
from typing import List, Optional

from backend.common.corpora_config import CorporaConfig
from backend.common.utils.aws import delete_many_from_s3
from backend.common.utils.result_notification import aws_batch_job_url_fmt_str, aws_sfn_url_fmt_str, notify_slack
from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetVersionId,
)
from backend.layers.processing import logger
from backend.portal.api.providers import get_business_logic

logger.configure_logging(level=logging.INFO)


def handle_failure(event: dict, context, delete_artifacts=True) -> None:
    logging.info(event)
    (
        dataset_version_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event(event)
    if "migrate" in execution_arn:
        logging.info(
            f"Skipping slack notification for {execution_arn} because failure is related to a schema migration"
        )
    else:
        trigger_slack_notification(
            dataset_version_id,
            collection_version_id,
            error_step_name,
            error_job_id,
            error_aws_regions,
            execution_arn,
        )
    update_dataset_processing_status_to_failed(dataset_version_id)
    if delete_artifacts:
        cleanup_artifacts(dataset_version_id, error_step_name)


def parse_event(event: dict):
    dataset_version_id = event.get("dataset_version_id")
    collection_version_id = event.get("collection_version_id")
    error_cause = event.get("error", {}).get("Cause", "")
    execution_arn = event.get("execution_id")
    try:
        error_cause_dict = json.loads(error_cause)
    except json.decoder.JSONDecodeError:
        error_step_name = None
        error_job_id = None
        error_aws_regions = None
    else:
        error_step_name = error_cause_dict.get("JobName")
        error_job_id = error_cause_dict.get("JobId")
        try:
            error_aws_regions = [
                i["Value"] for i in error_cause_dict["Container"]["Environment"] if i["Name"] == "AWS_DEFAULT_REGION"
            ][0]
        except KeyError:
            error_aws_regions = None
    return (
        dataset_version_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    )


def update_dataset_processing_status_to_failed(dataset_version_id: str) -> None:
    """
    This functions updates the processing status for a given dataset uuid to failed
    """
    with logger.LogSuppressed(Exception, message="Failed to update dataset processing status"):
        # If dataset not in db dont worry about updating its processing status
        get_business_logic().update_dataset_version_status(
            DatasetVersionId(dataset_version_id), DatasetStatusKey.PROCESSING, DatasetProcessingStatus.FAILURE
        )


def get_failure_slack_notification_message(
    dataset_version_id: Optional[str],
    collection_version_id: str,
    step_name: Optional[str],
    job_id: Optional[str],
    aws_region: str,
    execution_arn: str,
) -> dict:
    if dataset_version_id is None:
        logging.error("Dataset Version ID not found")
        dataset_version_id = "None"
        dataset = None
    else:
        dataset = get_business_logic().get_dataset_version(DatasetVersionId(dataset_version_id))
    if dataset is None:
        logging.error(f"Dataset version ID {dataset_version_id} not found")
        dataset_version_id = dataset_version_id + "(not found)"
        collection_owner, processing_status = "", ""
    else:
        collection_id = dataset.collection_id
        collection = get_business_logic().get_unpublished_collection_version_from_canonical(collection_id)
        if collection is None:
            logging.error(f"Collection {collection_id} not found")
            collection_owner = ""
        else:
            collection_owner = collection.owner
        processing_status = dataset.status.to_json(indent=2, sort_keys=True)
    batch_url = aws_batch_job_url_fmt_str.format(aws_region=aws_region, job_id=job_id)
    step_function_url = aws_sfn_url_fmt_str.format(aws_region=aws_region, execution_arn=execution_arn)
    collection_version_url = f"https://cellxgene.cziscience.com/collections/{collection_version_id}"

    # Generate canonical collection URL
    version = get_business_logic().get_collection_version(CollectionVersionId(collection_version_id))
    collection_url = f"https://cellxgene.cziscience.com/collections/{version.collection_id.id}"

    data = {
        "blocks": [
            {
                "type": "header",
                "text": {
                    "type": "plain_text",
                    "text": "Dataset failed to process:fire:",
                    "emoji": True,
                },
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! Please follow the triage steps: "
                    f"https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit"
                    f"#bookmark=id.3ofm47y0709y\n"
                    f"*Owner*: {collection_owner}\n"
                    f"*Collection URL*: {collection_url}\n"
                    f"*Collection Version URL*: {collection_version_url}\n"
                    f"*Batch Job ID*: <{batch_url}|{job_id}>\n"
                    f"*Step Function ARN*: <{step_function_url}|{execution_arn}>\n"
                    f"*Error Step*: {step_name}\n"
                    f"*Dataset Version ID*: {dataset_version_id}\n"
                    f"*Processing Status*:\n",
                },
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"```{processing_status}```",
                },
            },
        ]
    }
    return data


def trigger_slack_notification(
    dataset_version_id: Optional[str],
    collection_version_id: Optional[str],
    step_name: Optional[str],
    job_id: Optional[str],
    aws_region: str,
    execution_arn: str,
) -> None:
    with logger.LogSuppressed(Exception, message="Failed to send slack notification"):
        data = get_failure_slack_notification_message(
            dataset_version_id, collection_version_id, step_name, job_id, aws_region, execution_arn
        )
        # For these notifications, we should alert #single-cell-wrangling
        webhook = CorporaConfig().wrangling_slack_webhook
        notify_slack(data, webhook)


FAILED_ARTIFACT_CLEANUP_MESSAGE = "Failed to clean up artifacts."
FAILED_DATASET_CLEANUP_MESSAGE = "Failed to clean up datasets."
FAILED_CXG_CLEANUP_MESSAGE = "Failed to clean up cxgs."
FAILED_ATAC_DATASET_MESSAGE = "Failed to clean up ATAC datasets fragment files. artifact_id: {}"


def delete_atac_fragment_files(dataset_version_id: str) -> None:
    """Delete all atac fragment files and index files from S3 assocaited with the dataset version"""
    dataset = get_business_logic().get_dataset_version(DatasetVersionId(dataset_version_id))
    if not dataset:
        # If dataset not in db dont worry about deleting its files
        return

    if dataset.status.atac_status in [
        DatasetConversionStatus.COPIED,
        DatasetConversionStatus.SKIPPED,
        DatasetConversionStatus.NA,
    ]:
        # If the dataset is copied , we don't need to delete the files since they are part of another dataset.
        # If the dataset is skipped or NA, we don't need to delete the files since they are not created.
        return

    object_keys: List[str] = get_business_logic().get_atac_fragment_uris_from_dataset_version(dataset)
    for ok in object_keys:
        object_key = os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), ok)
        delete_and_catch_error("DATASETS_BUCKET", object_key, FAILED_ATAC_DATASET_MESSAGE.format(ok))


def delete_and_catch_error(bucket_name: str, object_key: str, error_message: str) -> None:
    with logger.LogSuppressed(Exception, message=error_message):
        bucket_name = os.environ[bucket_name]
        delete_many_from_s3(bucket_name, object_key)


def cleanup_artifacts(dataset_version_id: str, error_step_name: Optional[str] = None) -> None:
    """Clean up artifacts"""

    object_key = os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), dataset_version_id).strip("/")

    if error_step_name in ["validate_anndata", None]:
        delete_and_catch_error("ARTIFACT_BUCKET", object_key + "/", FAILED_ARTIFACT_CLEANUP_MESSAGE)
    if error_step_name in ["validate_atac", None]:
        delete_atac_fragment_files(dataset_version_id)

    delete_and_catch_error("DATASETS_BUCKET", object_key + ".", FAILED_DATASET_CLEANUP_MESSAGE)
    delete_and_catch_error("CELLXGENE_BUCKET", object_key + ".cxg/", FAILED_CXG_CLEANUP_MESSAGE)
