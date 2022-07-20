from flask import g, make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.api.v1.common import get_dataset_else_error


@dbconnect
def get(collection_id: str, curator_tag: str = None, dataset_id: str = None):
    db_session = g.db_session
    dataset = get_dataset_else_error(db_session, dataset_id, collection_id, curator_tag)
    status = dataset.processing_status.to_dict(
        remove_none=True, remove_attr=["dataset", "created_at", "updated_at", "id", "dataset_id"]
    )
    return make_response(jsonify(status), 200)
