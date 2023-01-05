import json
import logging
import os

from backend.common.utils.aws import delete_many_from_s3
from backend.common.utils.json import CustomJSONEncoder
from backend.common.utils.result_notification import format_failed_batch_issue_slack_alert, notify_slack
from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatusKey, DatasetVersionId
from backend.layers.persistence.persistence import DatabaseProvider

database_provider = DatabaseProvider()
business_logic = BusinessLogic(database_provider, None, None, None, None)


def handle_failure(event: dict, context) -> None:
    dataset_id = event["dataset_id"]
    trigger_slack_notification(dataset_id)
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
            DatasetVersionId(dataset_id), DatasetStatusKey.PROCESSING, DatasetProcessingStatus.FAILURE
        )
    # If dataset not in db dont worry about updating its processing status
    except Exception:
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


def trigger_slack_notification(dataset_id):
    dataset = business_logic.get_dataset_version(DatasetVersionId(dataset_id))
    if dataset is None:
        return
    collection_id = dataset.collection_id
    collection = business_logic.get_collection_version_from_canonical(collection_id)
    if collection is None:
        return
    collection_id, collection_owner = collection.version_id, collection.owner
    processing_status = dataset.status.to_json()

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
                    "text": f"Dataset processing job failed! @sc-oncall-eng\n"
                    f"*Owner*: {collection_owner}\n"
                    f"*Collection*: https://cellxgene.cziscience.com/collections/{collection_id}\n"
                    f"*Processing Status*:\n",
                },
            },
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"```{json.dumps(processing_status, cls=CustomJSONEncoder, indent=2, sort_keys=True)}```",
                },
            },
        ]
    }
    batch_alert_data = format_failed_batch_issue_slack_alert(data)
    logger = logging.getLogger(__name__)
    logger.info(batch_alert_data)
    notify_slack(data)
