from flask import make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.lambdas.api.v1.collection import post_collection_revision_common


@dbconnect
def post_collection_revision(collection_uuid: str, token_info: dict):
    collection_revision = post_collection_revision_common(collection_uuid, token_info)
    return make_response(jsonify({"revision_id": collection_revision.id}), 201)
