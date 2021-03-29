from flask import make_response, jsonify, g

from ....common.corpora_orm import CollectionVisibility
from ....common.entities.geneset import Geneset
from ....common.utils.exceptions import (
    ForbiddenHTTPException,
)


def delete(geneset_uuid: str, user: str):
    """
    Deletes an existing dataset or cancels an in progress upload.
    """
    db_session = g.db_session
    geneset = Geneset.get(db_session, geneset_uuid)
    accepted_response = "", 202
    if not geneset:
        return accepted_response
    if geneset.collection.owner != user:
        raise ForbiddenHTTPException()
    if geneset.collection_visibility == CollectionVisibility.PUBLIC:
        return make_response(jsonify("Can not delete a public geneset"), 405)
    geneset.delete()
    return accepted_response
