from flask import make_response

from ...db import get_db
from .....common.corpora_orm import CollectionVisibility
from .....common.entities import Collection
from .....common.utils.exceptions import ForbiddenHTTPException


def post(collection_uuid: str, user: str):
    session = get_db()
    collection = Collection.if_owner(session, collection_uuid, CollectionVisibility.PRIVATE, user)
    if not collection:
        raise ForbiddenHTTPException()
    collection.publish()
    return make_response({"collection_uuid": collection.id, "visibility": collection.visibility}, 202)
