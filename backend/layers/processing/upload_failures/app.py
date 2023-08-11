import json
import logging
import os
from typing import Optional

from backend.common.utils.aws import delete_many_from_s3
from backend.common.utils.result_notification import aws_batch_job_url_fmt_str, aws_sfn_url_fmt_str, notify_slack
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatusKey, DatasetVersionId
from backend.layers.processing import logger
from backend.portal.api.providers import get_business_logic

logger.configure_logging(level=logging.INFO)


def handle_failure(event: dict, context) -> None:
    logging.info(event)
    (
        dataset_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event(event)
    trigger_slack_notification(
        dataset_id, collection_version_id, error_step_name, error_job_id, error_aws_regions, execution_arn
    )
    update_dataset_processing_status_to_failed(dataset_id)
    cleanup_artifacts(dataset_id)


# write test cases using pytest to test the parse_event function


def parse_event(event: dict):
    dataset_id = event.get("dataset_id")
    collection_version_id = event.get("collection_id")
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
        dataset_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    )


def update_dataset_processing_status_to_failed(dataset_id: str) -> None:
    """
    This functions updates the processing status for a given dataset uuid to failed
    """
    with logger.LogSuppressed(Exception, message="Failed to update dataset processing status"):
        # If dataset not in db dont worry about updating its processing status
        get_business_logic().update_dataset_version_status(
            DatasetVersionId(dataset_id), DatasetStatusKey.PROCESSING, DatasetProcessingStatus.FAILURE
        )


def get_failure_slack_notification_message(
    dataset_id: Optional[str],
    collection_version_id: str,
    step_name: Optional[str],
    job_id: Optional[str],
    aws_region: str,
    execution_arn: str,
) -> dict:
    if dataset_id is None:
        logger.error("Dataset ID not found")  # type: ignore
        dataset_id = "None"
        dataset = None
    else:
        dataset = get_business_logic().get_dataset_version(DatasetVersionId(dataset_id))
    if dataset is None:
        logger.error(f"Dataset {dataset_id} not found")  # type: ignore
        dataset_id = dataset_id + "(not found)"
        collection_owner, processing_status = "", ""
    else:
        collection_id = dataset.collection_id
        collection = get_business_logic().get_unpublished_collection_version_from_canonical(collection_id)
        if collection is None:
            logger.error(f"Collection {collection_id} not found")  # type: ignore
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
                    f"*Dataset ID*: {dataset_id}\n"
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
    dataset_id: Optional[str],
    collection_version_id: Optional[str],
    step_name: Optional[str],
    job_id: Optional[str],
    aws_region: str,
    execution_arn: str,
) -> None:
    with logger.LogSuppressed(Exception, message="Failed to send slack notification"):
        data = get_failure_slack_notification_message(
            dataset_id, collection_version_id, step_name, job_id, aws_region, execution_arn  # type: ignore
        )
        notify_slack(data)


def cleanup_artifacts(dataset_id: str) -> None:
    """Clean up artifacts"""
    object_key = os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), dataset_id).strip("/")
    with logger.LogSuppressed(Exception, message="Failed to clean up artifacts."):
        artifact_bucket = os.environ["ARTIFACT_BUCKET"]
        delete_many_from_s3(artifact_bucket, object_key + "/")
    with logger.LogSuppressed(Exception, message="Failed to clean up datasets."):
        datasets_bucket = os.environ["DATASETS_BUCKET"]
        delete_many_from_s3(datasets_bucket, object_key + ".")
    with logger.LogSuppressed(Exception, message="Failed to clean up cxgs."):
        cellxgene_bucket = os.environ["CELLXGENE_BUCKET"]
        delete_many_from_s3(cellxgene_bucket, object_key + ".cxg/")
