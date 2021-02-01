from backend.corpora.common.entities import Dataset
from backend.corpora.common.corpora_orm import UploadStatus, DbDatasetProcessingStatus
from backend.corpora.common.utils.db_utils import processing_status_updater


def update_dataset_processing_status_to_failed(dataset_uuid, error=None) -> None:
    """
    This functions updates the processing status for a given dataset uuid to failed
    """
    try:
        status = {
            DbDatasetProcessingStatus.upload_progress: 0,
            DbDatasetProcessingStatus.upload_status: UploadStatus.FAILED,
            DbDatasetProcessingStatus.upload_message: str(error),
        }

        dataset = Dataset.get(dataset_uuid)
        processing_status_updater(dataset.processing_status.id, status)
    # If dataset not in db dont worry about updating its processing status
    except AttributeError:
        pass
