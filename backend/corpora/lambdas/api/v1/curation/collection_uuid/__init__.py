from flask import g

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.utils.http_exceptions import MethodNotAllowedException
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed
from backend.corpora.lambdas.api.v1.common import get_collection


@dbconnect
def delete_collection(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection(
        db_session,
        collection_uuid,
        visibility=CollectionVisibility.PRIVATE,
        owner=owner_or_allowed(token_info),
        include_tombstones=True,
    )
    if collection.revision_of:
        raise MethodNotAllowedException("Can only delete a private collection through curation API.")
    elif collection.tombstone:
        pass
    else:
        collection.delete()
    return "", 204
