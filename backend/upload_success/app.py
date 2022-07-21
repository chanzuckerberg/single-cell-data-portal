import logging

from backend.common.corpora_orm import ProcessingStatus, DbDatasetProcessingStatus
from backend.common.entities import Dataset
from backend.common.utils.db_helpers import processing_status_updater
from backend.common.utils.db_session import db_session_manager

logger = logging.getLogger(__name__)


def success_handler(event: dict, context) -> None:
    """
    Lambda function invoked by the ingestion step function that updates
    the processing status for the specified dataset to SUCCESS
    :param event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """
    dataset_id = event["dataset_id"]

    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id)
        success_status = {DbDatasetProcessingStatus.processing_status: ProcessingStatus.SUCCESS}
        processing_status_updater(session, dataset.processing_status.id, success_status)
