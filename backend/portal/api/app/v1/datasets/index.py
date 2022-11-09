from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.entities import Collection


@dbconnect
def get():
    db_session = g.db_session
    datasets = Collection.list_public_datasets_for_index(db_session)
    return make_response(jsonify(datasets), 200)
