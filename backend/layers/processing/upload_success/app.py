import logging

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatusKey, DatasetVersionId
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.upload_failures.app import handle_failure

database_provider = DatabaseProvider()
business_logic = BusinessLogic(database_provider, None, None, None, None)

logger = logging.getLogger("processing")


def success_handler(events: list, context) -> None:
    """
    Lambda function invoked by the ingestion step function that updates
    the processing status for the specified dataset to SUCCESS
    :param event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """
    cxg_event, seurat_event = events

    if "error" in cxg_event:
        handle_failure(cxg_event)
    else:
        if "error" in seurat_event:
            handle_failure(seurat_event)
        business_logic.update_dataset_version_status(
            DatasetVersionId(cxg_event["dataset_id"]), DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS
        )
