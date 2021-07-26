from backend.corpora.common.entities import Dataset
from backend.corpora.common.corpora_orm import UploadStatus, DbDatasetProcessingStatus, ProcessingStatus
from backend.corpora.common.utils.db_helpers import processing_status_updater
from backend.corpora.common.utils.db_session import db_session_manager


def update_dataset_processing_status_to_failed(dataset_uuid, error=None) -> None:
    """
    This functions updates the processing status for a given dataset uuid to failed
    """
    try:
        with db_session_manager() as session:
            dataset = Dataset.get(session, dataset_uuid)
            status = {
                DbDatasetProcessingStatus.processing_status: ProcessingStatus.FAILURE,
                DbDatasetProcessingStatus.upload_message: str(error),
            }
            processing_status_updater(session, dataset.processing_status.id, status)
    # If dataset not in db dont worry about updating its processing status
    except AttributeError:
        pass
