from flask import make_response, jsonify

from ....common.corpora_orm import DbDatasetProcessingStatus, UploadStatus
from ....common.entities import Dataset, Collection, DatasetAsset
from ....common.utils.db_utils import db_session, processing_status_updater
from ....common.utils.exceptions import (
    NotFoundHTTPException,
    ServerErrorHTTPException,
    ForbiddenHTTPException,
    MethodNotAllowedException,
)


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
def delete_dataset(dataset_uuid: str, user: str):
    """
    Deletes an existing dataset or cancels an in progress upload.
    """
    dataset = Dataset.get(dataset_uuid)
    if not dataset:
        raise ForbiddenHTTPException()
    if not Collection.if_owner(dataset.collection.id, dataset.collection.visibility, user):
        raise ForbiddenHTTPException()
    curr_status = dataset.processing_status
    if curr_status.upload_status is UploadStatus.UPLOADED:
        dataset.dataset_and_asset_deletion()
    elif curr_status.upload_status in [UploadStatus.UPLOADING, UploadStatus.WAITING]:
        status = {
            DbDatasetProcessingStatus.upload_progress: curr_status.upload_progress,
            DbDatasetProcessingStatus.upload_status: UploadStatus.CANCEL_PENDING,
        }
        processing_status_updater(dataset.processing_status.id, status)
        updated_status = Dataset.get(dataset_uuid).processing_status.to_dict()
        for remove in ["dataset", "created_at", "updated_at"]:
            updated_status.pop(remove)
        return make_response(jsonify(updated_status), 202)
    return "", 202
