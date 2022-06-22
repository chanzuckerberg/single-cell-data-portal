from flask import jsonify, g

from .common import get_access_type, reshape_for_curation_api
from .common import EntityColumns
from ......common.entities import Collection
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
    collections = Collection.list_collections_curation(g.db_session, EntityColumns.columns_for_collections, visibility)
    allowed_collections = []
    for collection in collections:
        access_type = get_access_type(collection, token_info)
        if access_type or visibility == CollectionVisibility.PUBLIC.name:
            reshape_for_curation_api(collection, access_type=access_type)
            allowed_collections.append(collection)

    return jsonify({"collections": allowed_collections})
