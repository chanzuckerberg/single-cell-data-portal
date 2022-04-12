import logging

from backend.corpora.common.corpora_orm import ProcessingStatus, DbDatasetProcessingStatus
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.utils.db_helpers import processing_status_updater
from backend.corpora.common.utils.db_session import db_session_manager

logger = logging.getLogger(__name__)


def success_handler(event: dict, context) -> None:
    """
    Lambda function invoked by the ingestion step function that updates
    the processing status for the specified dataset to SUCCESS
    :param event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """
    dataset_uuid = event["dataset_uuid"]

    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_uuid)
        success_status = {DbDatasetProcessingStatus.processing_status: ProcessingStatus.SUCCESS}
        processing_status_updater(session, dataset.processing_status.id, success_status)
