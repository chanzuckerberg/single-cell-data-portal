from flask import jsonify, g

from .common import reshape_for_curation_api_and_is_allowed, list_collections_curation
from .common import EntityColumns
from ......common.corpora_orm import CollectionVisibility
from ......common.utils.http_exceptions import UnauthorizedError
from backend.corpora.api_server.db import dbconnect


@dbconnect
def get(visibility: str, token_info: dict):
    """
    Collections index endpoint for Curation API. Only return Collection data for which the curator is authorized.
    :param visibility: the CollectionVisibility in string form
    :param token_info: access token info
    @return: Response
    """
    if not token_info and visibility == CollectionVisibility.PRIVATE.name:
        raise UnauthorizedError()
    collections = list_collections_curation(g.db_session, EntityColumns.columns_for_collections, visibility)
    allowed_collections = []
    for collection in collections:
        if reshape_for_curation_api_and_is_allowed(collection, token_info, preview=True):
            allowed_collections.append(collection)

    return jsonify({"collections": allowed_collections})
