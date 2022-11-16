from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.entities import Dataset
from backend.common.utils.http_exceptions import NotFoundHTTPException, ServerErrorHTTPException


@dbconnect
def post(dataset_id: str, asset_id: str):
    db_session = g.db_session
    # retrieve the dataset
    dataset = Dataset.get(db_session, dataset_id)
    if not dataset:
        raise NotFoundHTTPException(detail=f"'dataset/{dataset_id}' not found.")

    # retrieve the artifact
    asset = dataset.get_asset(asset_id)
    if not asset:
        raise NotFoundHTTPException(detail=f"'dataset/{dataset_id}/asset/{asset_id}' not found.")

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
            dataset_id=dataset_id,
            file_name=asset.filename,
            file_size=file_size,
            presigned_url=presigned_url,
        ),
        200,
    )
