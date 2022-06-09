from flask import g

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.utils.http_exceptions import MethodNotAllowedException
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed
from backend.corpora.lambdas.api.v1.common import get_collection


@dbconnect
def delete(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection(db_session, collection_uuid, owner=owner_or_allowed(token_info))
    if collection.visibility == CollectionVisibility.PUBLIC:
        raise MethodNotAllowedException("Cannot delete a public collection through API.")
    else:
        collection.delete()
    return "", 204
