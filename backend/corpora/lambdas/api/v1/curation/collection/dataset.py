import boto3
from flask import g, make_response

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_config import CorporaConfig
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.utils.exceptions import ForbiddenHTTPException, MethodNotAllowedException
from backend.corpora.lambdas.api.v1.collection import get_collection, _is_user_owner_or_allowed

sts_client = boto3.client("sts")

policy = """
{
    "Version": "2012-10-17",
    "Statement": [
        {
      	"Sid": "DataPortalUserUploadPolicy",
           	"Effect": "Allow",
           	"Action": [
           	        "s3:PutObject"
          	    ],
            "Resource": [
                    "arn:aws:s3:::{data_bucket}/{collection_id}/*"
                ]
        }
    ]
}
"""
duration = 3600


@dbconnect
def get_s3_credentials(collection_uuid, user):
    db_session = g.db_session
    config = CorporaConfig()
    # check if they own the collection.
    collection = get_collection(db_session, collection_uuid, include_tombstones=True)
    if not _is_user_owner_or_allowed(user, collection.owner):
        raise ForbiddenHTTPException()
    if collection.visibility != CollectionVisibility.PRIVATE.name:
        raise MethodNotAllowedException()
    credentials = sts_client.assume_role(
        RoleArn=config.curator_role_arn,
        RoleSessionName=user,
        Policy=policy.format(config.submission_bucket, collection_uuid),
        DurationSeconds=duration,
    )
    return make_response(credentials, 201)
