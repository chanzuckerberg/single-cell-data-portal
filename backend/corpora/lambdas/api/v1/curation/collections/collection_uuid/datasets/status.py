from flask import g, make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.lambdas.api.v1.common import get_dataset_else_invalid_parameter


@dbconnect
def get(collection_uuid: str, curator_tag: str = None, dataset_uuid: str = None):
    db_session = g.db_session
    dataset = get_dataset_else_invalid_parameter(db_session, dataset_uuid, collection_uuid, curator_tag)
    status = dataset.processing_status.to_dict(
        remove_none=True, remove_attr=["dataset", "created_at", "updated_at", "id", "dataset_id"]
    )
    return make_response(jsonify(status), 200)
