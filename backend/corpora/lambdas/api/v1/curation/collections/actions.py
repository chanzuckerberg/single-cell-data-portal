from flask import jsonify, g
from .common import add_collection_level_processing_status, reshape_for_curation_api
from .common import EntityColumns
from ...authorization import is_super_curator, owner_or_allowed
from ......common.corpora_orm import CollectionVisibility, DbCollection
from ......common.utils.http_exceptions import UnauthorizedError, ForbiddenHTTPException
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
    if visibility:
        filters.append(DbCollection.visibility == getattr(CollectionVisibility, visibility))
        if visibility == CollectionVisibility.PRIVATE.name:
            if not token_info:
                raise UnauthorizedError()
            else:
                result = owner_or_allowed(token_info)
                if result:
                    filters.append(DbCollection.owner == result)

    if curator:
        if not is_super_curator(token_info):
            raise ForbiddenHTTPException(detail="User is not authorized to use the curator parameter.")
        else:
            filters.append(DbCollection.curator_name == curator)

    db_session = g.db_session
    resp_collections = []
    for collection in db_session.query(DbCollection).filter(*filters).all():
        resp_collection = collection.to_dict_keep(EntityColumns.columns_for_collections)
        resp_collection["processing_status"] = add_collection_level_processing_status(collection)
        reshape_for_curation_api(db_session, resp_collection, token_info, preview=True)
        resp_collections.append(resp_collection)

    return jsonify({"collections": resp_collections})
