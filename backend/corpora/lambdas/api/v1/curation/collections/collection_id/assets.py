from flask import g, make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import DatasetAsset
from backend.corpora.common.utils.http_exceptions import NotFoundHTTPException
from backend.corpora.lambdas.api.v1.common import get_dataset_else_error


@dbconnect
def get(collection_id: str, dataset_id=None):
    db_session = g.db_session
    dataset = get_dataset_else_error(db_session, dataset_id, collection_id)

    # retrieve the artifact
    assets = dataset.get_assets()
    if not assets:
        raise NotFoundHTTPException(detail="No assets found. The dataset may still be processing.")

    asset_list = []
    error_flag = False
    for a in assets:
        asset = DatasetAsset(a)
        if asset.filetype == DatasetArtifactFileType.CXG:
            continue
        result = dict(filetype=asset.filetype, filename=asset.filename)

        # Retrieve S3 metadata
        filesize = asset.get_file_size()
        if not filesize:
            result["filesize"] = -1
            asset_list.append(result)
            error_flag = True
            continue
        else:
            result["filesize"] = filesize

        # Generate pre-signed URL
        presigned_url = asset.generate_file_url()
        if not presigned_url:
            result["presigned_url"] = "Not Found."
            error_flag = True
        else:
            result["presigned_url"] = presigned_url
        asset_list.append(result)

    response = dict(assets=asset_list)
    status_code = 202 if error_flag else 200
    return make_response(jsonify(**response), status_code)
