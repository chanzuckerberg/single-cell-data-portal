from flask import g, make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed
from backend.corpora.lambdas.api.v1.common import authorize_get_collection, validate_dataset_identifier


@dbconnect
def get(token_info: dict, collection_uuid: str, curator_tag: str = None, dataset_uuid: str = None):
    db_session = g.db_session
    dataset = validate_dataset_identifier(db_session, collection_uuid, dataset_uuid, curator_tag)
    authorize_get_collection(
        db_session,
        dataset.collection.id,
        owner=owner_or_allowed(token_info),
    )
    status = dataset.processing_status.to_dict(
        remove_none=True, remove_attr=["dataset", "created_at", "updated_at", "id", "dataset_id"]
    )
    return make_response(jsonify(status), 200)
