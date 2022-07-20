# Wraps tiledb_data.py to fit API expected responses

from flask import make_response, jsonify
from backend.corpora.lambdas.api.v1.collection_id.upload import upload_from_link
from backend.corpora.lambdas.api.v1.dataset import post_dataset_asset

from tiledb_data import TileDBData

location = "../../../../../tests/unit/backend/fixtures/test_tiledb/metadata"


def create_collection(body: dict, user: str):
    """/v1/collections POST"""
    db = TileDBData(location)
    from backend.corpora.lambdas.api.v1.collection import normalize_and_get_doi, get_publisher_metadata
    from backend.corpora.common.providers import crossref_provider

    errors = []
    doi = normalize_and_get_doi(body, errors)
    if errors:
        raise Exception(detail=errors)
    if doi is not None:
        provider = crossref_provider.CrossrefProvider()
        publisher_metadata = get_publisher_metadata(provider, doi)
    else:
        publisher_metadata = None
    id = db.create_collection(
        visibility="PRIVATE",
        name=body["name"],
        description=body["description"],
        owner=user,
        links=body.get("links", []),
        contact_name=body["contact_name"],
        contact_email=body["contact_email"],
        curator_name=body.get("curator_name", ""),
        publisher_metadata=publisher_metadata,
    )
    res = {"collection_id": id}
    return make_response(res, 200)


def list_published_collections(user_id: str, from_date: int, to_date: int):
    """/v1/collections GET"""
    db = TileDBData(location)
    colls = db.get_published_collections(user_id, from_date, to_date)
    res = {
        "collections": colls,
        "from_date": from_date,
        "to_date": to_date
    }
    return make_response(res, 200)


def list_published_collections_compact():
    """/v1/collections/index GET"""
    db = TileDBData(location)
    colls = db.get_published_collections()
    return make_response(jsonify(colls))


def delete_collection(id: str):
    """/v1/collections/{id} DELETE"""
    db = TileDBData(location)
    db.delete_collection(id)
    return make_response('', 204)


def get_collection(id: str):
    """/v1/collections/{id} GET"""
    db = TileDBData(location)
    coll = db.get_collection(id)
    datasets = db.get_datasets(coll)
    coll["datasets"] = datasets
    return make_response(jsonify(coll), 200)


def start_revision(id: str):
    """/v1/collections/{id} POST"""
    db = TileDBData(location)
    rev_id = db.create_revision(id)
    res = get_collection(rev_id)
    return make_response(res, 200)


def update_collection(id: str, body: dict):
    """/v1/collections/{id} PUT"""
    db = TileDBData(location)
    for key, val in body.items():
        db.edit_collection(id, key, val)
    res = get_collection(id)
    return make_response(res, 200)


def publish_collection(id: str):
    """/v1/collections/{id}/publish POST"""
    db = TileDBData(location)
    db.publish_collection(id)
    res = {
        "collection_id": id,
        "visibility": "PUBLIC"
    }
    return make_response(res, 202)


def upload_dataset(coll_id: str, body: dict, token_info: dict):
    """/v1/collections/{id}/upload-links POST"""
    # get the dataset data from the given url, manage and upload the artifact, etc.
    dataset_id = upload_from_link(coll_id, token_info, body["url"], curator_tag=body.get("curator_tag"))
    res = {"dataset_id": dataset_id}
    return make_response(res, 202)


def replace_dataset(coll_id: str,  body: dict, token_info: dict):
    """/v1/collections/{id}/upload-links PUT"""
    db = TileDBData(location)
    dataset_id = body.get("id")
    db.delete_dataset(coll_id, dataset_id)
    dataset_id = upload_from_link(
        coll_id,
        token_info,
        body.get("url", body.get("link")),
        dataset_id,
        curator_tag=body.get("curator_tag"),
    )
    res = {"dataset_id": dataset_id}
    return make_response(res, 202)


def list_published_datasets():
    """/v1/datasets/index GET"""
    db = TileDBData(location)
    datasets = db.get_published_datasets()
    return make_response(jsonify(datasets), 200)


def get_dataset_metadata(explorer_url: str):
    """/v1/datasets/meta GET"""
    # TODO
    return


def delete_dataset(coll_id: str, dataset_id: str):
    """/v1/datasets/{id} DELETE"""
    db = TileDBData(location)
    db.delete_dataset(coll_id, dataset_id)
    return make_response('', 202)


def request_asset_download(dataset_id: str, asset_id: str):
    """/v1/datasets/{dataset_id}/asset/{asset_id} POST"""
    res = post_dataset_asset(dataset_id, asset_id)
    return make_response(res, 200)


def list_dataset_assets(id: str):
    """/v1/datasets/{dataset_id}/assets GET"""
    db = TileDBData(location)
    dataset = db.get_dataset(id)
    assets = dataset["dataset_assets"]
    return make_response(jsonify(assets), 200)


def get_dataset_upload_status(id: str):
    """/v1/datasets/{dataset_id}/status GET"""
    db = TileDBData(location)
    dataset = db.get_dataset(id)
    status = dataset['processing_status']
    return make_response(status, 200)
