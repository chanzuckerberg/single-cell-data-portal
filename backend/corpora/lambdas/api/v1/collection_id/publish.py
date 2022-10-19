from flask import make_response, g

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Collection
from backend.corpora.common.utils.http_exceptions import ConflictException

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.utils.http_exceptions import ForbiddenHTTPException
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed

from backend.corpora.common.utils import cloudfront


@dbconnect
def post(collection_id: str, body: object, token_info: dict):
    db_session = g.db_session
    collection = Collection.get_collection(
        db_session,
        collection_id,
        CollectionVisibility.PRIVATE,
        owner=owner_or_allowed(token_info),
    )
    if not collection:
        raise ForbiddenHTTPException()
    if all([dataset.tombstone for dataset in collection.datasets]):
        raise ConflictException(detail="The collection must have a least one dataset.")
    collection_id = collection.revision_of if collection.revision_of else collection.id
    data_submission_policy_version = body["data_submission_policy_version"]
    collection.publish(data_submission_policy_version=data_submission_policy_version)
    cloudfront.create_invalidation_for_index_paths()
    return make_response({"collection_id": collection_id, "visibility": CollectionVisibility.PUBLIC}, 202)
