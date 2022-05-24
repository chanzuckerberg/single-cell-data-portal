import json

import sqlalchemy
from typing import Optional, List
from backend.corpora.common.providers import crossref_provider
from backend.corpora.common.providers.crossref_provider import CrossrefDOINotFoundException, CrossrefException

from flask import make_response, jsonify, g
from urllib.parse import urlparse
import re

import logging

from .common import get_collection
from ....common.corpora_orm import DbCollection, CollectionVisibility
from ....common.entities import Collection
from .authorization import is_user_owner_or_allowed, owner_or_allowed
from ....common.utils.http_exceptions import (
    InvalidParametersHTTPException,
    ConflictException,
)
from ....api_server.db import dbconnect


@dbconnect
def get_collections_list(from_date: int = None, to_date: int = None, token_info: Optional[dict] = None):
    db_session = g.db_session
    all_collections = Collection.list_attributes_in_time_range(
        db_session,
        from_date=from_date,
        to_date=to_date,
        list_attributes=[
            DbCollection.id,
            DbCollection.visibility,
            DbCollection.owner,
            DbCollection.created_at,
            DbCollection.revision_of,
        ],
    )

    collections = []
    for coll_dict in all_collections:
        visibility = coll_dict["visibility"]
        owner = coll_dict["owner"]
        if visibility == CollectionVisibility.PUBLIC:
            collections.append(dict(id=coll_dict["id"], created_at=coll_dict["created_at"], visibility=visibility.name))
        elif is_user_owner_or_allowed(token_info, owner):
            collections.append(
                dict(
                    id=coll_dict["id"],
                    created_at=coll_dict["created_at"],
                    visibility=visibility.name,
                    revision_of=coll_dict["revision_of"],
                )
            )

    result = {"collections": collections}
    if from_date:
        result["from_date"] = from_date
    if to_date:
        result["to_date"] = to_date

    return make_response(jsonify(result), 200)


@dbconnect
def get_collection_details(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection(db_session, collection_uuid, include_tombstones=True)
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
def get_collections_index():
    # TODO (ebezzi): this is very similar to `get_collections_list` above. Eventually they should be consolidated
    db_session = g.db_session

    filtered_collection = Collection.list_attributes_in_time_range(
        db_session,
        filters=[DbCollection.visibility == CollectionVisibility.PUBLIC],
        list_attributes=[
            DbCollection.id,
            DbCollection.name,
            DbCollection.published_at,
            DbCollection.revised_at,
            DbCollection.publisher_metadata,
        ],
    )

    # Remove entries where the value is None
    updated_collection = []
    for d in filtered_collection:
        updated_collection.append({k: v for k, v in d.items() if v is not None})

    return make_response(jsonify(updated_collection), 200)


@dbconnect
def post_collection_revision(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection(
        db_session,
        collection_uuid,
        visibility=CollectionVisibility.PUBLIC,
        owner=owner_or_allowed(token_info),
    )
    try:
        collection_revision = collection.create_revision()
    except sqlalchemy.exc.IntegrityError as ex:
        db_session.rollback()
        raise ConflictException() from ex
    result = collection_revision.reshape_for_api()

    result["access_type"] = "WRITE"
    return make_response(jsonify(result), 201)


doi_regex = re.compile(r"^10.\d{4,9}/[-._;()/:A-Z0-9]+$", flags=re.I)


def normalize_and_get_doi(body: dict, errors: list) -> Optional[str]:
    """
    1. Check for DOI uniqueness in the payload
    2. Normalizes it so that the DOI is always a link (starts with https://doi.org)
    3. Returns the newly normalized DOI
    """
    links = body.get("links", [])
    dois = [link for link in links if link["link_type"] == "DOI"]

    if not dois:
        return None

    # Verify that a single DOI exists
    if len(dois) > 1:
        errors.append({"link_type": "doi", "reason": "Can only specify a single DOI"})
        return None

    doi_node = dois[0]
    doi = doi_node["link_url"]

    parsed = urlparse(doi)
    if not parsed.scheme and not parsed.netloc:
        parsed_doi = parsed.path
        if not doi_regex.match(parsed_doi):
            errors.append({"link_type": "doi", "reason": "Invalid DOI"})
            return None
        doi_node["link_url"] = f"https://doi.org/{parsed_doi}"

    return doi


def get_publisher_metadata(provider, doi):
    """
    Retrieves publisher metadata from Crossref.
    """
    provider = crossref_provider.CrossrefProvider()
    try:
        return provider.fetch_metadata(doi)
    except CrossrefDOINotFoundException:
        # TODO: add an error message
        raise InvalidParametersHTTPException("DOI cannot be found on Crossref")
    except CrossrefException as e:
        logging.warning(f"CrossrefException on create_collection: {e}. Will ignore metadata.")
        return None


email_regex = re.compile(r"(.+)@(.+)\.(.+)")


def verify_collection_body(body: dict, errors: list, allow_none: bool = False) -> List[str]:
    result = email_regex.match(body.get("contact_email", ""))
    if not result and not allow_none:
        errors.append({"name": "contact_email", "reason": "Invalid format."})

    if not body.get("description") and not allow_none:  # Check if description is None or 0 length
        errors.append({"name": "description", "reason": "Cannot be blank."})

    if not body.get("name") and not allow_none:  # Check if name is None or 0 length
        errors.append({"name": "name", "reason": "Cannot be blank."})

    if not body.get("contact_name") and not allow_none:  # Check if contact_name is None or 0 length
        errors.append({"name": "contact_name", "reason": "Cannot be blank."})
    return errors


@dbconnect
def create_collection(body: dict, user: str):
    db_session = g.db_session
    errors = []
    verify_collection_body(body, errors)
    doi = normalize_and_get_doi(body, errors)
    if errors:
        raise InvalidParametersHTTPException(detail=json.dumps())
    if doi is not None:
        provider = crossref_provider.CrossrefProvider()
        publisher_metadata = get_publisher_metadata(provider, doi)
    else:
        publisher_metadata = None

    collection = Collection.create(
        db_session,
        visibility=CollectionVisibility.PRIVATE,
        name=body["name"],
        description=body["description"],
        owner=user,
        links=body.get("links", []),
        contact_name=body["contact_name"],
        contact_email=body["contact_email"],
        curator_name=body.get("curator_name", ""),
        publisher_metadata=publisher_metadata,
    )

    return make_response(jsonify({"collection_uuid": collection.id}), 201)


@dbconnect
def delete_collection(collection_uuid: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection(
        db_session,
        collection_uuid,
        owner=owner_or_allowed(token_info),
        include_tombstones=True,
    )
    if collection.visibility == CollectionVisibility.PUBLIC:
        revision = Collection.get_collection(
            db_session,
            revision_of=collection_uuid,
            owner=owner_or_allowed(token_info),
        )
        if revision:
            revision.delete()
        collection.tombstone_collection()
    else:
        collection.delete()
    return "", 204


@dbconnect
def update_collection(collection_uuid: str, body: dict, token_info: dict):
    db_session = g.db_session
    errors = []
    verify_collection_body(body, errors, allow_none=True)
    collection = get_collection(
        db_session,
        collection_uuid,
        visibility=CollectionVisibility.PRIVATE.name,
        owner=owner_or_allowed(token_info),
    )

    # Compute the diff between old and new DOI
    old_doi = collection.get_doi()
    new_doi = normalize_and_get_doi(body, errors)
    if errors:
        raise InvalidParametersHTTPException(detail=json.dumps())

    if old_doi and not new_doi:
        # If the DOI was deleted, remove the publisher_metadata field
        collection.update(publisher_metadata=None)
    elif new_doi != old_doi:
        # If the DOI has changed, fetch and update the metadata
        provider = crossref_provider.CrossrefProvider()
        publisher_metadata = get_publisher_metadata(provider, new_doi)
        body["publisher_metadata"] = publisher_metadata

    collection.update(**body)
    result = collection.reshape_for_api(tombstoned_datasets=True)
    result["access_type"] = "WRITE"
    return make_response(jsonify(result), 200)
