import boto3
import os
import sys

try:
    from backend.corpora.common.entities import Dataset
    from backend.corpora.common.utils.db_utils import db_session
    from backend.corpora.common.corpora_orm import UploadStatus, DbDatasetProcessingStatus
    from backend.corpora.dataset_processing.download import processing_status_updater
except ImportError:
    pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
    sys.path.insert(0, pkg_root)  # noqa
    from common.entities import Dataset
    from common.utils.db_utils import db_session
    from common.corpora_orm import UploadStatus, DbDatasetProcessingStatus
    from dataset_processing.download import processing_status_updater


def delete_many_from_s3(bucket_name, dataset_uuid):
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    bucket.objects.filter(Prefix=f"{dataset_uuid}/").delete()


@db_session
def update_dataset_processing_status_to_failed(dataset_uuid, error=None):
    try:

        status = {
                DbDatasetProcessingStatus.upload_progress: 0,
                DbDatasetProcessingStatus.upload_status: UploadStatus.FAILED,
                DbDatasetProcessingStatus.upload_message: str(error)
            }
        dataset = Dataset.get(dataset_uuid)
        processing_status_updater(dataset.processing_status.id, status)
    # If dataset not in db dont worry about updating its processing status
    except AttributeError:
        pass

