import logging

from backend.corpora.common.corpora_orm import ProcessingStatus, DbDatasetProcessingStatus
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_helpers import processing_status_updater
from backend.corpora.common.utils.db_session import db_session_manager

logger = logging.getLogger(__name__)


def dataset_submission_handler(event: dict, context) -> None:
    """
    Lambda function invoked when a dataset is uploaded to the cellxgene-dataset-submissions-{dev,staging,prod} S3 bucket
    :param event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """

    # TODO: implement me