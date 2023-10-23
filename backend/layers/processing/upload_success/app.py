import logging

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatusKey, DatasetVersionId
from backend.layers.persistence.persistence import DatabaseProvider

database_provider = DatabaseProvider()
business_logic = BusinessLogic(database_provider, None, None, None, None)

logger = logging.getLogger("processing")


def success_handler(event: dict, context) -> None:
    """
    Lambda function invoked by the ingestion step function that updates
    the processing status for the specified dataset to SUCCESS
    :param event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """
    dataset_id = event["dataset_id"]

    business_logic.update_dataset_version_status(
        DatasetVersionId(dataset_id), DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS
    )
