from flask import make_response, jsonify, g

from ....common.corpora_orm import DbDatasetProcessingStatus, UploadStatus
from ....common.entities import Dataset, Collection
from ....common.utils.db_session import processing_status_updater
from ....common.utils.exceptions import (
    NotFoundHTTPException,
    ServerErrorHTTPException,
    ForbiddenHTTPException,
    MethodNotAllowedException,
)


def post_dataset_asset(dataset_uuid: str, asset_uuid: str):
    db_session = g.db_session
    # retrieve the dataset
    dataset = Dataset.get(db_session, dataset_uuid)
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


def get_status(dataset_uuid: str, user: str):
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_uuid)
    if not Collection.if_owner(db_session, dataset.collection.id, dataset.collection.visibility, user):
        raise ForbiddenHTTPException()
    status = dataset.processing_status.to_dict(remove_none=True)
    for remove in ["dataset", "created_at", "updated_at"]:
        status.pop(remove)
    return make_response(jsonify(status), 200)


def delete_dataset(dataset_uuid: str, user: str):
    """
    Cancels an inprogress upload.
    """
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_uuid)
    if not dataset:
        raise ForbiddenHTTPException()
    if not Collection.if_owner(db_session, dataset.collection.id, dataset.collection.visibility, user):
        raise ForbiddenHTTPException()
    curr_status = dataset.processing_status
    if curr_status.upload_status is UploadStatus.UPLOADED:
        raise MethodNotAllowedException(f"'dataset/{dataset_uuid}' upload is complete and can not be cancelled.")
    status = {
        DbDatasetProcessingStatus.upload_progress: curr_status.upload_progress,
        DbDatasetProcessingStatus.upload_status: UploadStatus.CANCEL_PENDING,
    }
    processing_status_updater(db_session, dataset.processing_status.id, status)
    db_session.refresh(dataset.db_object)
    updated_status = Dataset.get(db_session, dataset_uuid).processing_status.to_dict()
    for remove in ["dataset", "created_at", "updated_at"]:
        updated_status.pop(remove)
    return make_response(jsonify(updated_status), 202)
