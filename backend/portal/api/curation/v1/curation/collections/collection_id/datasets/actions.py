from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.corpora_orm import DbCollection, CollectionVisibility, ProcessingStatus
from backend.common.entities import Dataset
from backend.common.utils.http_exceptions import (
    NotFoundHTTPException,
    MethodNotAllowedException,
)
from backend.portal.api.app.v1.authorization import owner_or_allowed
from backend.portal.api.app.v1.collection_id.upload_links import upload_from_link
from backend.portal.api.collections_common import (
    get_dataset_else_error,
    delete_dataset_common,
    get_collection_else_forbidden,
)
from backend.portal.api.curation.v1.curation.collections.common import EntityColumns, reshape_dataset_for_curation_api


@dbconnect
def get(collection_id: str, dataset_id: str = None):
    db_session = g.db_session
    if not db_session.query(DbCollection.id).filter(DbCollection.id == collection_id).first():
        raise NotFoundHTTPException("Collection not found!")
    dataset = get_dataset_else_error(db_session, dataset_id, collection_id)
    response_body = reshape_dataset_for_curation_api(dataset.to_dict_keep(EntityColumns.columns_for_dataset))
    return make_response(jsonify(response_body), 200)


@dbconnect
def delete(token_info: dict, collection_id: str, dataset_id: str = None):
    db_session = g.db_session
    dataset = get_dataset_else_error(db_session, dataset_id, collection_id, include_tombstones=True)
    delete_dataset_common(db_session, dataset, token_info)
    return "", 202


@dbconnect
def post(token_info: dict, collection_id: str):
    db_session = g.db_session
    collection = get_collection_else_forbidden(db_session, collection_id, owner=owner_or_allowed(token_info))
    if collection.visibility != CollectionVisibility.PRIVATE:
        raise MethodNotAllowedException("Collection must be PRIVATE Collection, or a revision of a PUBLIC Collection.")
    dataset = Dataset.create(
        db_session, collection=collection, processing_status={"processing_status": ProcessingStatus.INITIALIZED}
    )
    return make_response(jsonify({"id": dataset.id}), 201)


def put(collection_id: str, dataset_id: str, body: dict, token_info: dict):
    upload_from_link(
        collection_id,
        token_info,
        body.get("url", body.get("link")),
        dataset_id,
    )
    return "", 202
