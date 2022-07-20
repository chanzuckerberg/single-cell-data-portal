import json
import logging

import boto3
from flask import g, request, make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.api.v1.authorization import owner_or_allowed, is_super_curator
from backend.corpora.api.v1.common import get_collection_else_forbidden, get_dataset_else_error
from backend.corpora.api.v1.dataset import delete_dataset_common


@dbconnect
def delete_dataset(token_info: dict, collection_id: str, curator_tag: str = None, dataset_id: str = None):
    db_session = g.db_session
    dataset = get_dataset_else_error(db_session, dataset_id, collection_id, curator_tag, include_tombstones=True)
    delete_dataset_common(db_session, dataset, token_info)
    return "", 202


sts_client = boto3.client("sts")
logger = logging.getLogger(__name__)
duration = 43200


@dbconnect
def post_s3_credentials(collection_id: str, token_info: dict):
    db_session = g.db_session
    config = CorporaConfig()
    # Raise an error if they are not allowed to modify the collection.
    get_collection_else_forbidden(
        db_session, collection_id, visibility=CollectionVisibility.PRIVATE.name, owner=owner_or_allowed(token_info)
    )
    user_id = token_info["sub"]
    if is_super_curator(token_info):
        upload_key_prefix = f"super/{collection_id}/"
    else:
        upload_key_prefix = f"{user_id}/{collection_id}/"
    parameters = dict(
        RoleArn=config.curator_role_arn,
        RoleSessionName=user_id.replace("|", "-"),
        WebIdentityToken=request.headers["Authorization"].split(" ")[1],
        DurationSeconds=duration,
    )
    logger.info(json.dumps(parameters))
    response = sts_client.assume_role_with_web_identity(**parameters)
    response["UploadKeyPrefix"] = f"{upload_key_prefix}"
    response["Bucket"] = config.submission_bucket
    return make_response(jsonify(response), 200)
