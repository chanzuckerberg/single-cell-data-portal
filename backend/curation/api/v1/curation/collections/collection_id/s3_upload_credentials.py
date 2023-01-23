import json
import logging

import boto3
from flask import jsonify, make_response, request

from backend.common.corpora_config import CorporaConfig
from backend.common.utils.http_exceptions import ForbiddenHTTPException
from backend.layers.auth.user_info import UserInfo
from backend.curation.api.v1.curation.collections.common import (
    get_infered_collection_version_else_forbidden,
    is_owner_or_allowed_else_forbidden,
)

sts_client = boto3.client("sts")

logger = logging.getLogger(__name__)
duration = 43200


def get(collection_id: str, token_info: dict):
    config = CorporaConfig()
    user_info = UserInfo(token_info)
    # TODO: Since this method only works on private collections, I think we should just accept the version_id here
    collection_version = get_infered_collection_version_else_forbidden(collection_id)
    is_owner_or_allowed_else_forbidden(collection_version, user_info)
    if collection_version.published_at:
        raise ForbiddenHTTPException()
    if user_info.is_super_curator():
        upload_key_prefix = f"super/{collection_id}/"
    else:
        upload_key_prefix = f"{user_info.user_id}/{collection_id}/"
    parameters = dict(
        RoleArn=config.curator_role_arn,
        RoleSessionName=user_info.user_id.replace("|", "-"),
        WebIdentityToken=request.headers["Authorization"].split(" ")[1],
        DurationSeconds=duration,
    )
    logger.info(json.dumps(parameters))
    response = sts_client.assume_role_with_web_identity(**parameters)
    response["UploadKeyPrefix"] = f"{upload_key_prefix}"
    response["Bucket"] = config.submission_bucket
    return make_response(jsonify(response), 200)
