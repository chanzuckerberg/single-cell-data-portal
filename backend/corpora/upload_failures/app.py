import os

from backend.corpora.common.corpora_orm import DbDatasetProcessingStatus, ProcessingStatus
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.aws import delete_many_from_s3
from backend.corpora.common.utils.db_helpers import processing_status_updater
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.slack import dataset_processing_slack_notification


def handle_failure(event: dict, context) -> None:
    dataset_id = event["dataset_id"]
    dataset_processing_slack_notification(dataset_id)
    object_key = os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), dataset_id).strip("/")
    delete_many_from_s3(os.environ["ARTIFACT_BUCKET"], object_key)
    cellxgene_bucket = f"hosted-cellxgene-{os.environ['DEPLOYMENT_STAGE']}"
    if os.getenv("CELLXGENE_BUCKET"):
        cellxgene_bucket = os.getenv("CELLXGENE_BUCKET")
    delete_many_from_s3(cellxgene_bucket, object_key)
    update_dataset_processing_status_to_failed(dataset_id, event["error"]["Cause"])


def update_dataset_processing_status_to_failed(dataset_id, error=None) -> None:
    """
    This functions updates the processing status for a given dataset uuid to failed
    """
    try:
        with db_session_manager() as session:
            dataset = Dataset.get(session, dataset_id)
            status = {
                DbDatasetProcessingStatus.processing_status: ProcessingStatus.FAILURE,
                DbDatasetProcessingStatus.upload_message: str(error),
            }
            processing_status_updater(session, dataset.processing_status.id, status)
    # If dataset not in db dont worry about updating its processing status
    except AttributeError:
        pass
