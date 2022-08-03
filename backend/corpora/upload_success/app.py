import logging

from backend.corpora.common.corpora_orm import ProcessingStatus
from backend.corpora.common.entities import db
from backend.corpora.dataset_processing.process import update_db

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

    dataset = db.get_dataset(dataset_id)
    success_status = {"processing_status": ProcessingStatus.SUCCESS.name}
    update_db(dataset_id, processing_status=success_status)
