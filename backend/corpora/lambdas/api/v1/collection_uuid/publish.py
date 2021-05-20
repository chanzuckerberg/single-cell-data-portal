from flask import make_response, g

from .....common.corpora_orm import CollectionVisibility
from .....common.entities import Collection
from .....api_server.db import dbconnect
from .....common.utils.exceptions import ForbiddenHTTPException


@dbconnect
def post(collection_uuid: str, user: str):
    db_session = g.db_session
    collection = Collection.get_collection(db_session, collection_uuid, CollectionVisibility.PRIVATE, owner=user)
    if not collection:
        raise ForbiddenHTTPException()
    collection.publish()
    return make_response({"collection_uuid": collection.id, "visibility": collection.visibility}, 202)
