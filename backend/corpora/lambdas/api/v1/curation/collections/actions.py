from flask import jsonify, g
from .common import reshape_for_curation_api_and_is_allowed, add_collection_level_processing_status
from .common import EntityColumns
from ...authorization import is_super_curator
from ......common.corpora_orm import CollectionVisibility, DbCollection
from ......common.utils.http_exceptions import UnauthorizedError
from backend.corpora.api_server.db import dbconnect


@dbconnect
def get(visibility: str, token_info: dict, curator: str = None):
    """
    Collections index endpoint for Curation API. Only return Collection data for which the curator is authorized.
    :param visibility: the CollectionVisibility in string form
    :param token_info: access token info
    @return: Response
    """
    if visibility == CollectionVisibility.PRIVATE.name:
        if any((not token_info, not is_super_curator(token_info) and curator not in ("self", None))):
            raise UnauthorizedError()

    filters = [DbCollection.tombstone == False]  # noqa
    if visibility == CollectionVisibility.PUBLIC.name:
        filters.append(DbCollection.visibility == CollectionVisibility.PUBLIC)
    elif visibility == CollectionVisibility.PRIVATE.name:
        filters.append(DbCollection.visibility == CollectionVisibility.PRIVATE)

    if token_info and curator == "self":  # curator wants own collections
        filters.append(DbCollection.owner == token_info["sub"])
    elif curator:  # user want collections from a specific curator
        filters.append(DbCollection.curator_name == curator)

    db_session = g.db_session
    resp_collections = []
    for collection in db_session.query(DbCollection).filter(*filters).all():
        resp_collection = collection.to_dict_keep(EntityColumns.columns_for_collections)
        resp_collection["processing_status"] = add_collection_level_processing_status(collection)
        if reshape_for_curation_api_and_is_allowed(db_session, resp_collection, token_info, preview=True):
            resp_collections.append(resp_collection)

    return jsonify({"collections": resp_collections})
