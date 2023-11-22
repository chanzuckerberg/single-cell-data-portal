import contextlib
import json
import logging
import os
from typing import Optional

from backend.common.corpora_config import CorporaConfig
from backend.common.utils.aws import delete_many_from_s3
from backend.common.utils.result_notification import aws_batch_job_url_fmt_str, aws_sfn_url_fmt_str, notify_slack
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatusKey, DatasetVersionId
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
    trigger_slack_notification(
        dataset_version_id, collection_version_id, error_step_name, error_job_id, error_aws_regions, execution_arn
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
                    "text": f"Dataset processing job failed! @sc-oncall-eng please follow the [triage steps](https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit#bookmark=id.3ofm47y0709y)\n"
                    f"*Owner*: {collection_owner}\n"
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
    with contextlib.suppress(Exception):
        data = get_failure_slack_notification_message(
            dataset_version_id, collection_version_id, step_name, job_id, aws_region, execution_arn
        )
        # For these notifications, we should alert #single-cell-wrangling
        webhook = CorporaConfig().wrangling_slack_webhook
        notify_slack(data, webhook)


FAILED_ARTIFACT_CLEANUP_MESSAGE = "Failed to clean up artifacts."
FAILED_DATASET_CLEANUP_MESSAGE = "Failed to clean up datasets."
FAILED_CXG_CLEANUP_MESSAGE = "Failed to clean up cxgs."


def cleanup_artifacts(dataset_version_id: str, error_step_name: Optional[str] = None) -> None:
    """Clean up artifacts"""
    object_key = os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), dataset_version_id).strip("/")
    if not error_step_name or error_step_name in ["validate", "download"]:
        with logger.LogSuppressed(Exception, message=FAILED_ARTIFACT_CLEANUP_MESSAGE):
            artifact_bucket = os.environ["ARTIFACT_BUCKET"]
            delete_many_from_s3(artifact_bucket, object_key + "/")
    with logger.LogSuppressed(Exception, message=FAILED_DATASET_CLEANUP_MESSAGE):
        datasets_bucket = os.environ["DATASETS_BUCKET"]
        delete_many_from_s3(datasets_bucket, object_key + ".")
    with logger.LogSuppressed(Exception, message=FAILED_CXG_CLEANUP_MESSAGE):
        cellxgene_bucket = os.environ["CELLXGENE_BUCKET"]
        delete_many_from_s3(cellxgene_bucket, object_key + ".cxg/")
