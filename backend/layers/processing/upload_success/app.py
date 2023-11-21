import logging

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatusKey, DatasetVersionId
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.processing.upload_failures.app import handle_failure

database_provider = DatabaseProvider()
business_logic = BusinessLogic(database_provider, None, None, None, None)

logger = logging.getLogger("processing")


def success_handler(events: dict, context) -> None:
    """
    Lambda function invoked by the ingestion step function that updates
    the processing status for the specified dataset to SUCCESS
    :param event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """
    cxg_job, seurat_job = events["cxg_job"], events["seurat_job"]
    cxg_job["execution_id"], seurat_job["execution_id"] = events["execution_id"], events["execution_id"]

    if cxg_job.get("error"):
        handle_failure(cxg_job, context)
        exit(1)
    elif seurat_job.get("error"):
        handle_failure(seurat_job, context, delete_artifacts=False)
        exit(1)
    else:
        business_logic.update_dataset_version_status(
            DatasetVersionId(cxg_job["dataset_version_id"]),
            DatasetStatusKey.PROCESSING,
            DatasetProcessingStatus.SUCCESS,
        )
