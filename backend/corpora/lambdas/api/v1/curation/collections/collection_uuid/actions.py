from flask import g, jsonify

from ..common import EntityColumns, get_access_type, reshape_for_curation_api
from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Collection
from backend.corpora.common.utils.http_exceptions import MethodNotAllowedException, NotFoundHTTPException
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed
from backend.corpora.lambdas.api.v1.common import get_collection_else_forbidden


@dbconnect
def delete(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection_else_forbidden(db_session, collection_uuid, owner=owner_or_allowed(token_info))
    if collection.visibility == CollectionVisibility.PUBLIC:
        raise MethodNotAllowedException("Cannot delete a public collection through API.")
    else:
        collection.delete()
    return "", 204


@dbconnect
def get(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    collection = Collection.get_collection(db_session, collection_uuid, include_tombstones=True)
    if not collection:
        raise NotFoundHTTPException
    collection_response: dict = collection.to_dict_keep(EntityColumns.columns_for_collection_uuid)
    access_type = get_access_type(collection_response, token_info, uuid_provided=True)
    reshape_for_curation_api(collection_response, access_type=access_type)

    return jsonify(collection_response)
