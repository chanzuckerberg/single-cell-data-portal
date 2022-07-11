# Wraps tiledb_data.py to fit API expected responses

from flask import make_response, jsonify

from tiledb_data import TileDBData

location = "/Users/ragarwal/code/single-cell-data-portal/tests/unit/backend/fixtures/test_tiledb/metadata"

def create_collection(
    name: str = "",
    description: str = "",
    owner: str = "",
    contact_name: str = "",
    contact_email: str = "",
    curator_name: str = "",
    links: list = None,
):
    db = TileDBData(location)
    id = db.create_collection(name, description, owner, contact_name, contact_email, curator_name, links)
    res = {"collection_id": id}
    return make_response(jsonify(res), 200)

def list_published_collections(user_id, from_date, to_date):
    db = TileDBData(location)
    colls = db.get_published_collections(user_id, from_date, to_date)
    res = {
        "collections": colls,
        "from_date": from_date,
        "to_date": to_date
    }
    return make_response(jsonify(res), 200)

def list_published_collections_compact():
    db = TileDBData(location)
    colls = db.get_all_collections()
    return make_response(jsonify(colls))

