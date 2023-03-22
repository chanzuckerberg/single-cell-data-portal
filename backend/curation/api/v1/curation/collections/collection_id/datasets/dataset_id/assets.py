from flask import jsonify, make_response

from backend.common.utils.http_exceptions import NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import (
    allowed_dataset_asset_types,
    get_inferred_collection_version_else_forbidden,
    get_inferred_dataset_version,
)
from backend.portal.api.providers import get_business_logic


def get(collection_id: str, dataset_id=None):
    business_logic = get_business_logic()

    dataset = get_inferred_dataset_version(dataset_id)

    if dataset is None:
        get_inferred_collection_version_else_forbidden(collection_id)
        raise NotFoundHTTPException(detail="Dataset not found.")

    assets = dataset.artifacts
    if not assets:
        raise NotFoundHTTPException(detail="No assets found. The dataset may still be processing.")

    asset_list = []
    error_flag = False
    for asset in assets:
        if asset.type not in allowed_dataset_asset_types:
            continue
        download_data = business_logic.get_dataset_artifact_download_data(dataset.version_id, asset.id)
        if download_data.file_size is None:
            error_flag = True
            download_data.file_size = -1
        if not download_data.presigned_url:
            error_flag = True
            download_data.presigned_url = "Not Found."

        result = {
            "filename": download_data.file_name,
            "filesize": download_data.file_size,
            "filetype": download_data.file_type.upper(),
            "presigned_url": download_data.presigned_url,
        }

        asset_list.append(result)

    response = asset_list
    status_code = 202 if error_flag else 200
    return make_response(jsonify(response), status_code)
