import json
import logging

import boto3
from flask import g, request, make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.http_exceptions import InvalidParametersHTTPException
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed
from backend.corpora.lambdas.api.v1.common import get_collection
from backend.corpora.lambdas.api.v1.dataset import delete_dataset_common


@dbconnect
def delete_dataset(token_info: dict, collection_uuid: str, curator_tag: str = None, dataset_uuid=None):
    db_session = g.db_session
    if dataset_uuid:
        dataset = Dataset.get(db_session, dataset_uuid, include_tombstones=True)
    elif curator_tag:
        dataset = Dataset.get_dataset_from_curator_tag(db_session, collection_uuid, curator_tag)
    else:
        raise InvalidParametersHTTPException()
    delete_dataset_common(db_session, dataset, token_info)
    return "", 202


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
