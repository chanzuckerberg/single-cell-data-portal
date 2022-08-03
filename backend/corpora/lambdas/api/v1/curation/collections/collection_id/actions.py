from flask import g, jsonify

from backend.corpora.lambdas.api.v1.collection import update_collection
from ..common import EntityColumns, add_collection_level_processing_status
from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Collection
from backend.corpora.common.utils.http_exceptions import MethodNotAllowedException, NotFoundHTTPException
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed
from backend.corpora.lambdas.api.v1.curation.collections.common import reshape_for_curation_api_and_is_allowed
from backend.corpora.lambdas.api.v1.common import get_collection_else_forbidden


@dbconnect
def delete(collection_id: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection_else_forbidden(db_session, collection_id, owner=owner_or_allowed(token_info))
    if collection.visibility == CollectionVisibility.PUBLIC:
        raise MethodNotAllowedException(detail="Cannot delete a public collection through API.")
    else:
        collection.delete()
    return "", 204


@dbconnect
def get(collection_id: str, token_info: dict):
    db_session = g.db_session
    collection = Collection.get_collection(db_session, collection_id, include_tombstones=False)
    if not collection:
        raise NotFoundHTTPException
    collection_response: dict = collection.to_dict_keep(EntityColumns.columns_for_collection_id)
    collection_response["processing_status"] = add_collection_level_processing_status(collection)
    reshape_for_curation_api_and_is_allowed(collection_response, token_info, id_provided=True)

    return jsonify(collection_response)


def patch(collection_id: str, body: dict, token_info: dict):
    return update_collection(collection_id, body, token_info)
