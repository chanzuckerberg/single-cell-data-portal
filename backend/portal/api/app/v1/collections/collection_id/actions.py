from flask import jsonify, make_response, g

from backend.api_server.db import dbconnect
from backend.common.corpora_orm import CollectionVisibility
from backend.common.entities import Collection
from backend.common.utils.http_exceptions import InvalidParametersHTTPException
from backend.portal.api.app.v1.authorization import is_user_owner_or_allowed, owner_or_allowed
from backend.portal.api.app.v1.collections.actions import get_doi_link_node
from backend.portal.api.app.v1.common import portal_get_normalized_doi_url
from backend.portal.api.collections_common import (
    get_collection_else_forbidden,
    post_collection_revision_common,
    get_collection_and_verify_body,
)
from backend.portal.api.collections_common import get_publisher_metadata


@dbconnect
def delete(collection_id: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection_else_forbidden(db_session, collection_id, owner=owner_or_allowed(token_info))
    if collection.visibility == CollectionVisibility.PUBLIC:
        revision = Collection.get_collection(
            db_session,
            revision_of=collection_id,
            owner=owner_or_allowed(token_info),
        )
        if revision:
            revision.delete()
        collection.tombstone_collection()
    else:
        collection.delete()
    return "", 204


@dbconnect
def get(collection_id: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection_else_forbidden(db_session, collection_id, include_tombstones=True)
    if collection.tombstone:
        result = ""
        response = 410
    else:
        get_tombstone_datasets = (
            is_user_owner_or_allowed(token_info, collection.owner)
            and collection.visibility == CollectionVisibility.PRIVATE
        )
        result = collection.reshape_for_api(get_tombstone_datasets)
        response = 200
        result["access_type"] = "WRITE" if is_user_owner_or_allowed(token_info, collection.owner) else "READ"
    return make_response(jsonify(result), response)


@dbconnect
def post(collection_id: str, token_info: dict):
    collection_revision = post_collection_revision_common(collection_id, token_info)
    result = collection_revision.reshape_for_api()
    result["access_type"] = "WRITE"
    return make_response(jsonify(result), 201)


@dbconnect
def put(collection_id: str, body: dict, token_info: dict):
    db_session = g.db_session
    collection, errors = get_collection_and_verify_body(db_session, collection_id, body, token_info)
    # Compute the diff between old and new DOI
    old_doi = collection.get_doi()
    new_doi_url = None
    if new_doi_node := get_doi_link_node(body, errors):
        if new_doi_url := portal_get_normalized_doi_url(new_doi_node, errors):
            new_doi_node["link_url"] = new_doi_url
    if old_doi and not new_doi_url:
        # If the DOI was deleted, remove the publisher_metadata field
        collection.update(publisher_metadata=None)
    elif new_doi_url != old_doi:
        # If the DOI has changed, fetch and update the metadata
        publisher_metadata = get_publisher_metadata(new_doi_url, errors)
        body["publisher_metadata"] = publisher_metadata
    if errors:
        raise InvalidParametersHTTPException(detail=errors)
    collection.update(**body)
    result = collection.reshape_for_api(tombstoned_datasets=True)
    result["access_type"] = "WRITE"
    return make_response(jsonify(result), 200)
