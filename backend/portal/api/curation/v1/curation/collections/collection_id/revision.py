from flask import make_response, jsonify

from backend.api_server.db import dbconnect
from backend.portal.api.collections_common import post_collection_revision_common


@dbconnect
def post(collection_id: str, token_info: dict):
    collection_revision = post_collection_revision_common(collection_id, token_info)
    return make_response(jsonify({"id": collection_revision.id}), 201)
