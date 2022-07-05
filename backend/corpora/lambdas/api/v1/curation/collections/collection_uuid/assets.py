from flask import g, make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import DatasetAsset
from backend.corpora.common.utils.http_exceptions import NotFoundHTTPException
from backend.corpora.lambdas.api.v1.common import get_dataset_else_error


@dbconnect
def get(collection_uuid: str, curator_tag: str = None, dataset_uuid=None):
    db_session = g.db_session
    dataset = get_dataset_else_error(db_session, dataset_uuid, collection_uuid, curator_tag)

    # retrieve the artifact
    assets = dataset.get_assets()
    if not assets:
        raise NotFoundHTTPException("No assets found. The dataset may still be processing.")

    asset_list = []
    error_flag = False
    for a in assets:
        asset = DatasetAsset(a)
        if asset.filetype == DatasetArtifactFileType.CXG:
            continue
        result = dict(file_type=asset.filetype, file_name=asset.filename)

        # Retrieve S3 metadata
        file_size = asset.get_file_size()
        if not file_size:
            result["file_size"] = -1
            asset_list.append(result)
            error_flag = True
            continue
        else:
            result["file_size"] = file_size

        # Generate pre-signed URL
        presigned_url = asset.generate_file_url()
        if not presigned_url:
            result["presigned_url"] = "Not Found."
            error_flag = True
        else:
            result["presigned_url"] = presigned_url
        asset_list.append(result)

    response = dict(dataset_uuid=dataset.id, assets=asset_list)
    if dataset.curator_tag:
        response["curator_tag"] = dataset.curator_tag
    status_code = 202 if error_flag else 200
    return make_response(jsonify(**response), status_code)
