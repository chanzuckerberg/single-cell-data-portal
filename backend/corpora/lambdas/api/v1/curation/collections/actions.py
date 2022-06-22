from flask import jsonify, g

from .common import reshape_for_curation_api_and_is_allowed
from ......common.entities import Collection
from ......common.corpora_orm import CollectionVisibility
from ......common.utils.http_exceptions import UnauthorizedError
from backend.corpora.api_server.db import dbconnect


@dbconnect
def get_collections(visibility: str, token_info: dict):
    """
    Collections index endpoint for Curation API. Only return Collection data for which the curator is authorized.
    @param visibility: the CollectionVisibility in string form
    @param token_info: access token info
    @return: Response
    """
    if not token_info and visibility == CollectionVisibility.PRIVATE.name:
        raise UnauthorizedError()
    collections = Collection.list_collections_for_curation(g.db_session, visibility)
    allowed_collections = []
    for collection in collections:
        if reshape_for_curation_api_and_is_allowed(collection, token_info):
            allowed_collections.append(collection)

    return jsonify({"collections": allowed_collections})
