from flask import make_response

from .....common.corpora_orm import CollectionVisibility
from .....common.utils.db_utils import db_session_manager
from .....common.entities import Collection
from .....common.utils.exceptions import ForbiddenHTTPException


def post(collection_uuid: str, user: str):
    with db_session_manager() as session:
        collection = Collection.if_owner(session, collection_uuid, CollectionVisibility.PRIVATE, user)
        if not collection:
            raise ForbiddenHTTPException()
        collection.publish()
        return make_response({"collection_uuid": collection.id, "visibility": collection.visibility}, 202)
