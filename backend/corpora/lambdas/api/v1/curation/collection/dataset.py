import json
import boto3
import logging
from flask import g, make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.lambdas.api.v1.common import get_collection
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed

sts_client = boto3.client("sts")
logger = logging.getLogger(__name__)


def create_policy(data_bucket: str, collection_id: str) -> str:
    policy = {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Sid": "DataPortalUserUploadPolicy",
                "Effect": "Allow",
                "Action": ["s3:PutObject"],
                "Resource": [f"arn:aws:s3:::{data_bucket}/{collection_id}/*"],
            }
        ],
    }
    return json.dumps(policy)


duration = 3600


@dbconnect
def post_s3_credentials(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    config = CorporaConfig()
    # Raise an error if they are not allowed to modify the collection.
    get_collection(
        db_session, collection_uuid, visibility=CollectionVisibility.PRIVATE.name, owner=owner_or_allowed(token_info)
    )

    parameters = dict(
        RoleArn=config.curator_role_arn,
        RoleSessionName=token_info["sub"].replace("|", "-"),
        Policy=create_policy(config.submission_bucket, collection_uuid),
        DurationSeconds=duration,
    )
    logger.info(json.dumps(parameters))
    credentials = sts_client.assume_role(**parameters)
    return make_response(jsonify(credentials), 200)
