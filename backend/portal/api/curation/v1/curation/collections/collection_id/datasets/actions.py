from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.corpora_orm import CollectionVisibility, ProcessingStatus
from backend.common.entities import Dataset
from backend.common.utils.http_exceptions import MethodNotAllowedException
from backend.portal.api.app.v1.authorization import owner_or_allowed
from backend.portal.api.collections_common import get_collection_else_forbidden


@dbconnect
def post(token_info: dict, collection_id: str):
    db_session = g.db_session
    collection = get_collection_else_forbidden(db_session, collection_id, owner=owner_or_allowed(token_info))
    if collection.visibility != CollectionVisibility.PRIVATE:
        raise MethodNotAllowedException("Collection must be PRIVATE Collection, or a revision of a PUBLIC Collection.")
    dataset = Dataset.create(
        db_session, collection=collection, processing_status={"processing_status": ProcessingStatus.INITIALIZED}
    )
    return make_response(jsonify({"id": dataset.id}), 201)
