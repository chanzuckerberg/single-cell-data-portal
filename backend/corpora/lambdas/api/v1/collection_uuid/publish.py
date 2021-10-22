from flask import make_response, g

from .....common.corpora_orm import CollectionVisibility
from .....common.entities import Collection
from .....common.utils.exceptions import ConflictException

from .....api_server.db import dbconnect
from .....common.utils.exceptions import ForbiddenHTTPException
from backend.corpora.lambdas.api.v1.collection import _owner_or_allowed


@dbconnect
def post(collection_uuid: str, body: object, user: str):
    db_session = g.db_session
    collection = Collection.get_collection(
        db_session,
        collection_uuid,
        CollectionVisibility.PRIVATE,
        owner=_owner_or_allowed(user),
    )
    if not collection:
        raise ForbiddenHTTPException()
    if all([dataset.tombstone for dataset in collection.datasets]):
        raise ConflictException(detail="The collection must have a least one dataset.")

    data_submission_policy_version = body["data_submission_policy_version"]
    if data_submission_policy_version is None:
        raise ConflictException(detail="data_submission_policy_version not found")

    collection.publish(data_submission_policy_version=data_submission_policy_version)
    return make_response({"collection_uuid": collection.id, "visibility": collection.visibility}, 202)
