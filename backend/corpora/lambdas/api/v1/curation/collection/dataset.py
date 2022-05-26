import json
import boto3
import logging
from flask import g, make_response, jsonify, request

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.lambdas.api.v1.common import get_collection
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed

sts_client = boto3.client("sts")
logger = logging.getLogger(__name__)

duration = 43200


@dbconnect
def post_s3_credentials(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    config = CorporaConfig()
    # Raise an error if they are not allowed to modify the collection.
    get_collection(
        db_session, collection_uuid, visibility=CollectionVisibility.PRIVATE.name, owner=owner_or_allowed(token_info)
    )
    user_id = token_info["sub"]
    upload_key_prefix = f"{user_id}/{collection_uuid}"
    parameters = dict(
        RoleArn=config.curator_role_arn,
        RoleSessionName=user_id.replace("|", "-"),
        WebIdentityToken=request.headers["Authorization"].split(" ")[1],
        DurationSeconds=duration,
    )
    logger.info(json.dumps(parameters))
    response = sts_client.assume_role_with_web_identity(**parameters)
    response["UploadKeyPrefix"] = f"{upload_key_prefix}/"
    response["Bucket"] = config.submission_bucket
    return make_response(jsonify(response), 200)
