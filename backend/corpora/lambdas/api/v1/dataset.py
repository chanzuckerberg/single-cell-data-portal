from flask import make_response, jsonify

from ....common.entities import Dataset
from ....common.utils.db_utils import db_session
from ....common.utils.exceptions import NotFoundHTTPException, ServerErrorHTTPException


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
