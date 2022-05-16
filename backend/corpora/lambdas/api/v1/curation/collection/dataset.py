import json
import boto3
import logging
from flask import g, make_response, jsonify, request

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.utils.corpora_constants import CorporaConstants
from backend.corpora.lambdas.api.v1.common import get_collection
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed

sts_client = boto3.client("sts")
logger = logging.getLogger(__name__)


def create_policy(bucket: str, upload_path: str) -> str:
    policy = {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Sid": "DataPortalUserUploadPolicy",
                "Effect": "Allow",
                "Action": ["s3:PutObject"],
                "Resource": [f"arn:aws:s3:::{bucket}/{upload_path}/*"],
            }
        ],
    }
    return json.dumps(policy)


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
    if CorporaConstants.SUPER_CURATOR_SCOPE in token_info["scope"]:
        upload_path = f"{CorporaConstants.SUPER_CURATOR_NAME}/{collection_uuid}"
    else:
        upload_path = f"{user_id}/{collection_uuid}"
    parameters = dict(
        RoleArn=config.curator_role_arn,
        RoleSessionName=user_id.replace("|", "-"),
        Policy=create_policy(config.submission_bucket, upload_path),
        WebIdentityToken=request.headers["Authorization"][6:],
        DurationSeconds=duration,
    )
    logger.info(json.dumps(parameters))
    response = sts_client.assume_role_with_web_identity(**parameters)
    response["UploadPath"] = f"{upload_path}/"
    response["Bucket"] = config.submission_bucket
    return make_response(jsonify(response), 200)
