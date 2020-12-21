from flask import make_response

from .....common.corpora_orm import CollectionVisibility
from .....common.utils.db_utils import db_session
from .....common.entities import Collection
from .....common.utils.exceptions import ForbiddenHTTPException


@db_session()
def post(collection_uuid: str, user: str):
    private_collection = Collection.if_owner(collection_uuid, CollectionVisibility.PRIVATE, user)
    if not private_collection:
        raise ForbiddenHTTPException
    public_collection = Collection.publish(private_collection)
    return make_response({"collection_uuid": public_collection.id, "visibility": public_collection.visibility}, 202)
