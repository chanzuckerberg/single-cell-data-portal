from flask import make_response, jsonify, g

from backend.common.corpora_orm import CollectionVisibility
from backend.common.entities.geneset import Geneset
from backend.api_server.db import dbconnect_and_ddtrace
from backend.common.utils.authorization_checks import is_user_owner_or_allowed
from backend.common.utils.http_exceptions import (
    ForbiddenHTTPException,
)


@dbconnect_and_ddtrace
def delete(geneset_id: str, token_info: dict):
    """
    Deletes an existing geneset
    """
    db_session = g.db_session
    geneset = Geneset.get(db_session, geneset_id)
    accepted_response = "", 202
    if not geneset:
        return accepted_response
    if not is_user_owner_or_allowed(token_info["sub"], token_info.get("scope"), geneset.collection.owner):
        raise ForbiddenHTTPException()
    if geneset.collection.visibility == CollectionVisibility.PUBLIC:
        return make_response(jsonify("Cannot delete a public geneset"), 405)
    geneset.delete()
    return accepted_response
