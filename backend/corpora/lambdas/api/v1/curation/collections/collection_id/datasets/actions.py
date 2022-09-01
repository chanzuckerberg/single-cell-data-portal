from flask import g, make_response, jsonify

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_orm import CollectionVisibility, DbCollection
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.http_exceptions import (
    InvalidParametersHTTPException,
    ConflictException,
    NotFoundHTTPException,
)
from backend.corpora.common.utils.regex import validate_curator_tag
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed
from backend.corpora.lambdas.api.v1.common import (
    get_dataset_else_error,
    get_collection_else_forbidden,
    delete_dataset_common,
)
from backend.corpora.lambdas.api.v1.curation.collections.common import EntityColumns, reshape_dataset_for_curation_api


@dbconnect
def patch(token_info: dict, collection_id: str, body: dict, curator_tag: str = None, dataset_id: str = None):
    db_session = g.db_session
    dataset = get_dataset_else_error(db_session, dataset_id, collection_id, curator_tag)
    get_collection_else_forbidden(
        db_session, collection_id, visibility=CollectionVisibility.PRIVATE, owner=owner_or_allowed(token_info)
    )
    tag = body["curator_tag"]
    try:
        validate_curator_tag(tag)
    except ValueError as ex:
        raise InvalidParametersHTTPException(detail=ex.args[0])

    # Check if the curator_tag is unique across datasets in the collection.
    if dataset.curator_tag != tag:
        if conflict := Dataset.get(db_session, collection_id=collection_id, curator_tag=tag):
            raise ConflictException(
                detail=f"Curator_tags must be unique within a collection. Dataset={conflict.id} is using "
                f"curator_tag={tag}"
            )
        else:
            dataset.update(curator_tag=tag)
    return make_response("", 204)


@dbconnect
def get(collection_id: str, curator_tag: str = None, dataset_id: str = None):
    db_session = g.db_session
    if not db_session.query(DbCollection.id).filter(DbCollection.id == collection_id).first():
        raise NotFoundHTTPException("Collection not found!")
    dataset = get_dataset_else_error(db_session, dataset_id, collection_id, curator_tag)
    response_body = reshape_dataset_for_curation_api(dataset.to_dict_keep(EntityColumns.columns_for_dataset))
    return make_response(jsonify(response_body), 200)


@dbconnect
def delete(token_info: dict, collection_id: str, curator_tag: str = None, dataset_id: str = None):
    db_session = g.db_session
    dataset = get_dataset_else_error(db_session, dataset_id, collection_id, curator_tag, include_tombstones=True)
    delete_dataset_common(db_session, dataset, token_info)
    return "", 202
