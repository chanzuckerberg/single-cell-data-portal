from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.entities import Dataset
from backend.common.utils.http_exceptions import ForbiddenHTTPException
from backend.portal.api.app.v1.authorization import owner_or_allowed
from backend.portal.api.collections_common import get_collection_else_forbidden


@dbconnect
def get(dataset_id: str, token_info: dict):
    db_session = g.db_session
    dataset = Dataset.get(db_session, dataset_id)
    if not dataset:
        raise ForbiddenHTTPException()
    get_collection_else_forbidden(
        db_session,
        dataset.collection.id,
        owner=owner_or_allowed(token_info),
    )
    status = dataset.processing_status.to_dict(remove_none=True, remove_attr=["dataset", "created_at", "updated_at"])
    return make_response(jsonify(status), 200)
