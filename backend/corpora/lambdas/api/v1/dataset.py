from flask import make_response, jsonify

from ....common.corpora_orm import DbDatasetProcessingStatus, UploadStatus
from ....common.entities import Dataset, Collection
from ....common.utils.db_utils import db_session
from ....common.utils.exceptions import NotFoundHTTPException, ServerErrorHTTPException, ForbiddenHTTPException
from ....dataset_processing.download import processing_status_updater


@db_session()
def post_dataset_asset(dataset_uuid: str, asset_uuid: str):

    # retrieve the dataset
    dataset = Dataset.get(dataset_uuid)
    if not dataset:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}' not found.")

    # retrieve the artifact
    asset = dataset.get_asset(asset_uuid)
    if not asset:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}/asset/{asset_uuid}' not found.")

    # Retrieve S3 metadata
    file_size = asset.get_file_size()
    if not file_size:
        raise ServerErrorHTTPException()

    # Generate pre-signed URL
    presigned_url = asset.generate_file_url()
    if not presigned_url:
        raise ServerErrorHTTPException()

    return make_response(
        jsonify(
            dataset_id=dataset_uuid,
            file_name=asset.filename,
            file_size=file_size,
            presigned_url=presigned_url,
        ),
        200,
    )


@db_session()
def get_status(dataset_uuid: str, user: str):
    dataset = Dataset.get(dataset_uuid)
    if not Collection.if_owner(dataset.collection.id, dataset.collection.visibility, user):
        raise ForbiddenHTTPException()
    status = dataset.processing_status.to_dict(remove_none=True)
    for remove in ["dataset", "created_at", "updated_at"]:
        status.pop(remove)
    return make_response(jsonify(status), 200)


@db_session()
def cancel_dataset_download(dataset_uuid: str, user: str):
    dataset = Dataset.get(dataset_uuid)
    if not dataset:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}' not found.")
    if not Collection.if_owner(dataset.collection.id, dataset.collection.visibility, user):
        raise ForbiddenHTTPException()
    status = dataset.processing_status.to_dict()
    if status.upload_status is UploadStatus.UPLOADED:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}' upload is complete and can not be cancelled.")
    status = {
        DbDatasetProcessingStatus.upload_progress: status.upload_progress,
        DbDatasetProcessingStatus.upload_status: UploadStatus.CANCEL_PENDING,
    }
    processing_status_updater(dataset.processing_status.id, status)
    return make_response(jsonify(status), 202)
