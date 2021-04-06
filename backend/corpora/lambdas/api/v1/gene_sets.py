from flask import make_response, jsonify, g

from ....common.corpora_orm import CollectionVisibility
from ....common.entities.geneset import Geneset
from ....common.utils.db_session import dbconnect
from ....common.utils.exceptions import (
    ForbiddenHTTPException,
)


@dbconnect
def delete(geneset_uuid: str, user: str):
    """
    Deletes an existing geneset
    """
    db_session = g.db_session
    geneset = Geneset.get(db_session, geneset_uuid)
    accepted_response = "", 202
    if not geneset:
        return accepted_response
    if geneset.collection.owner != user:
        raise ForbiddenHTTPException()
    if geneset.collection_visibility == CollectionVisibility.PUBLIC:
        return make_response(jsonify("Cannot delete a public geneset"), 405)
    geneset.delete()
    return accepted_response
