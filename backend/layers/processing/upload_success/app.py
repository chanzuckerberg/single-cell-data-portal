import logging

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatusKey, DatasetVersionId
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.upload_failures.app import (
    handle_failure,
    parse_event,
    trigger_slack_notification,
    update_dataset_processing_status_to_failed,
)

database_provider = DatabaseProvider()
business_logic = BusinessLogic(database_provider, None, None, None, None)

logger = logging.getLogger("processing")


def success_handler(events: dict, context) -> None:
    """
    Lambda function invoked by the ingestion step function that updates
    the processing status for the specified dataset to SUCCESS
    :param event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """
    cxg_job, seurat_job = events["cxg_job"], events["seurat_job"]
    cxg_job["execution_id"], seurat_job["execution_id"] = events["execution_id"], events["execution_id"]

    if cxg_job.get("error"):
        handle_failure(cxg_job, context)
    elif seurat_job.get("error"):
        # Same as handle_failure except do not delete artifacts (because cxg succeeded)
        (
            dataset_id,
            collection_version_id,
            error_step_name,
            error_job_id,
            error_aws_regions,
            error_cause,
            execution_arn,
        ) = parse_event(seurat_job)
        update_dataset_processing_status_to_failed(dataset_id)
        trigger_slack_notification(
            dataset_id, collection_version_id, error_step_name, error_job_id, error_aws_regions, execution_arn
        )
    else:
        business_logic.update_dataset_version_status(
            DatasetVersionId(cxg_job["dataset_id"]), DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS
        )
