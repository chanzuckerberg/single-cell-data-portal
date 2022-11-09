from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.entities import Dataset


@dbconnect
def get(dataset_id: str):
    db_session = g.db_session
    # retrieve the dataset
    dataset = Dataset.get(db_session, dataset_id)
    assets = dataset.get_assets()
    return make_response(jsonify(assets=assets))
