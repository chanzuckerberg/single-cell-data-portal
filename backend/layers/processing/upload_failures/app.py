import json
import os

from backend.common.utils.aws import delete_many_from_s3
from backend.common.utils.result_notification import dataset_processing_slack_notification
from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatusKey, DatasetVersionId

database_provider = DatabaseProvider()
business_logic = BusinessLogic(database_provider, None, None, None, None)


def handle_failure(event: dict, context) -> None:
    dataset_id = event["dataset_id"]
    dataset_processing_slack_notification(dataset_id)
    update_dataset_processing_status_to_failed(dataset_id, event["error"]["Cause"])

    error_step_name = get_error_step_name(event)
    object_key = os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), dataset_id).strip("/")
    # clean up artifacts depending on error; default to full clean-up if error step is unknown
    if not error_step_name or error_step_name == "download-validate":
        delete_many_from_s3(os.environ["ARTIFACT_BUCKET"], object_key)
    if not error_step_name or error_step_name == "cxg":
        cellxgene_bucket = os.getenv("CELLXGENE_BUCKET", default=f"hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}")
        delete_many_from_s3(cellxgene_bucket, object_key)


def update_dataset_processing_status_to_failed(dataset_id, error=None) -> None:
    """
    This functions updates the processing status for a given dataset uuid to failed
    """
    try:
        business_logic.update_dataset_version_status(
            DatasetVersionId(dataset_id), 
            DatasetStatusKey.PROCESSING, 
            DatasetProcessingStatus.FAILURE
        )
    # If dataset not in db dont worry about updating its processing status
    except AttributeError:
        pass


def get_error_step_name(event: dict) -> str:
    """
    This function derives the name of the step that failed
    """
    error_cause_dict = json.loads(event["error"]["Cause"])
    error_job_name = None
    if "JobName" in error_cause_dict:
        error_job_name = error_cause_dict["JobName"]
    return error_job_name
