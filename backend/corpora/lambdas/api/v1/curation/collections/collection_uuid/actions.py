from flask import g, jsonify

from ..common import EntityColumns
from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Collection
from backend.corpora.common.utils.http_exceptions import MethodNotAllowedException, NotFoundHTTPException
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed
from backend.corpora.lambdas.api.v1.common import get_collection
from backend.corpora.lambdas.api.v1.curation.collections.common import reshape_for_curation_api_and_is_allowed


@dbconnect
def delete(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection(db_session, collection_uuid, owner=owner_or_allowed(token_info))
    if collection.visibility == CollectionVisibility.PUBLIC:
        raise MethodNotAllowedException("Cannot delete a public collection through API.")
    else:
        collection.delete()
    return "", 204


@dbconnect
def get_collection_uuid(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    collection = Collection.get_collection(db_session, collection_uuid, include_tombstones=True)
    if not collection:
        raise NotFoundHTTPException
    collection_response: dict = collection.to_dict_keep(EntityColumns.columns_for_collection_uuid)

    reshape_for_curation_api_and_is_allowed(collection_response, token_info, allow_access=True)

    return jsonify(collection_response)
