from flask import jsonify, g
from .common import reshape_for_curation_api
from ...authorization import is_super_curator, owner_or_allowed
from ...collection import create_collection
from ......common.corpora_orm import CollectionVisibility, DbCollection
from ......common.utils.http_exceptions import ForbiddenHTTPException, InvalidParametersHTTPException
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
                raise ForbiddenHTTPException(detail="Not authorized to query for PRIVATE collection.")
            else:
                owner_filter = owner_or_allowed(token_info)
                if owner_filter:  # None means the user is a super curator and don't need to filter by owner.
                    filters.append(DbCollection.owner == owner_filter)

    if curator:
        if not is_super_curator(token_info):
            raise ForbiddenHTTPException(detail="Not authorized to use the curator query parameter.")
        else:
            filters.append(DbCollection.curator_name == curator)

    db_session = g.db_session
    resp_collections = []
    for collection in db_session.query(DbCollection).filter(*filters).all():
        resp_collection = reshape_for_curation_api(db_session, collection, token_info, preview=True)
        resp_collections.append(resp_collection)

    return jsonify(resp_collections)


def post(body: dict, user: str):
    try:
        return create_collection(body, user)
    except InvalidParametersHTTPException as ex:
        ex.ext = dict(invalid_parameters=ex.detail)
        ex.detail = InvalidParametersHTTPException._default_detail
        raise ex
