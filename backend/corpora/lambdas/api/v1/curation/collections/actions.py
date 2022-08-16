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
    :param curator: the name of the curator to filter the collections by.

    @return: Response
    """
    filters = [DbCollection.tombstone == False]  # noqa

    if visibility == CollectionVisibility.PRIVATE.name and not token_info:
        raise UnauthorizedError()
    elif visibility:
        filters.append(DbCollection.visibility == getattr(CollectionVisibility, visibility))

    if curator and not is_super_curator(token_info):
        raise UnauthorizedError()
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
