import logging

from backend.corpora.common.corpora_orm import ProcessingStatus
from backend.corpora.dataset_processing.process import update_db

logger = logging.getLogger(__name__)


def success_handler(event: dict, context) -> None:
    """
    Lambda function invoked by the ingestion step function that updates the processing status
    for the specified dataset to SUCCESS
    :param event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """
    dataset_uuid = event["dataset_uuid"]

    try:
        update_db(dataset_uuid, processing_status=dict(processing_status=ProcessingStatus.SUCCESS))
    except Exception:
        logger.exception(f"Dataset {dataset_uuid}: Failed to update processing_status to SUCCESS")
