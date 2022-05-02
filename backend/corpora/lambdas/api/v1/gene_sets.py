from flask import make_response, jsonify, g

from ....common.corpora_orm import CollectionVisibility
from ....common.entities.geneset import Geneset
from ....api_server.db import dbconnect
from ....common.utils.authorization_checks import is_user_owner_or_allowed
from ....common.utils.exceptions import (
    ForbiddenHTTPException,
)


@dbconnect
def delete(geneset_uuid: str, token_info: dict):
    """
    Deletes an existing geneset
    """
    db_session = g.db_session
    geneset = Geneset.get(db_session, geneset_uuid)
    accepted_response = "", 202
    if not geneset:
        return accepted_response
    if not is_user_owner_or_allowed(token_info["sub"], geneset.collection.owner):
        raise ForbiddenHTTPException()
    if geneset.collection.visibility == CollectionVisibility.PUBLIC:
        return make_response(jsonify("Cannot delete a public geneset"), 405)
    geneset.delete()
    return accepted_response
