import sqlalchemy
from sqlalchemy.orm import Session
from typing import Optional
from backend.common.providers.crossref_provider import (
    CrossrefDOINotFoundException,
    CrossrefException,
    CrossrefProvider,
)

from flask import make_response, jsonify, g
from urllib.parse import urlparse
import re

import logging

from backend.portal.api.app.v1.common import get_collection_else_forbidden
from backend.common.corpora_orm import DbCollection, CollectionVisibility, ProjectLinkType
from backend.common.entities import Collection
from backend.portal.api.app.v1.authorization import is_user_owner_or_allowed, owner_or_allowed
from backend.common.utils.http_exceptions import (
    InvalidParametersHTTPException,
    ForbiddenHTTPException,
)
from backend.api_server.db import dbconnect
from backend.common.utils.regex import CONTROL_CHARS, DOI_REGEX_COMPILED, CURIE_REFERENCE_REGEX


@dbconnect
def get(from_date: int = None, to_date: int = None, token_info: Optional[dict] = None):
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
def get_collections_index():
    # TODO (ebezzi): this is very similar to `get` above. Eventually they should be consolidated
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


def post_collection_revision_common(collection_id: str, token_info: dict):
    db_session = g.db_session
    collection = get_collection_else_forbidden(
        db_session,
        collection_id,
        visibility=CollectionVisibility.PUBLIC,
        owner=owner_or_allowed(token_info),
    )
    try:
        collection_revision = collection.create_revision()
    except sqlalchemy.exc.IntegrityError:
        db_session.rollback()
        raise ForbiddenHTTPException("A revision is already in progess.")
    return collection_revision





def get_publisher_metadata(doi: str, errors: list) -> Optional[dict]:
    """
    Retrieves publisher metadata from Crossref.
    """
    provider = CrossrefProvider()
    try:
        return provider.fetch_metadata(doi)
    except CrossrefDOINotFoundException:
        errors.append({"link_type": ProjectLinkType.DOI.name, "reason": "DOI cannot be found on Crossref"})
    except CrossrefException as e:
        logging.warning(f"CrossrefException on create_collection: {e}. Will ignore metadata.")
        return None


email_regex = re.compile(r"(.+)@(.+)\.(.+)")


def verify_collection_links(body: dict, errors: list) -> None:
    def _error_message(i: int, _url: str) -> dict:
        return {"name": f"links[{i}]", "reason": "Invalid URL.", "value": _url}

    for index, link in enumerate(body.get("links", [])):
        if link["link_type"] == ProjectLinkType.DOI.name:
            continue
        url = link["link_url"]
        try:
            result = urlparse(url.strip())
        except ValueError:
            errors.append(_error_message(index, url))
        if not all([result.scheme, result.netloc]):
            errors.append(_error_message(index, url))


control_char_re = re.compile(CONTROL_CHARS)


def verify_collection_body(body: dict, errors: list) -> None:
    def check(key) -> bool:
        if key in body.keys():
            if not body[key]:
                errors.append({"name": key, "reason": "Cannot be blank."})
            elif key == "name" and control_char_re.search(body[key]):
                errors.append({"name": key, "reason": "Invalid characters detected."})
            else:
                return body[key]

    contact_email = check("contact_email")
    if contact_email:
        result = email_regex.match(contact_email)
        if not result:
            errors.append({"name": "contact_email", "reason": "Invalid format."})

    check("description")
    check("name")
    check("contact_name")

    verify_collection_links(body, errors)


def get_doi_link_node(body: dict, errors: list) -> Optional[dict]:
    links = body.get("links", [])
    dois = [link for link in links if link["link_type"] == ProjectLinkType.DOI.name]

    if not dois:
        return None

    # Verify that a single DOI exists
    if len(dois) > 1:
        errors.append({"link_type": ProjectLinkType.DOI.name, "reason": "Can only specify a single DOI"})
        return None

    return dois[0]


def curation_get_normalized_doi_url(doi: str, errors: list) -> Optional[str]:
    # Regex below is adapted from https://bioregistry.io/registry/doi 'Pattern for CURIES'
    if not re.match(CURIE_REFERENCE_REGEX, doi):
        errors.append({"name": ProjectLinkType.DOI.value, "reason": "DOI must be a CURIE reference."})
        return None
    return f"https://doi.org/{doi}"


def portal_get_normalized_doi_url(doi_node: dict, errors: list) -> Optional[str]:
    """
    1. Check for DOI uniqueness in the payload
    2. Normalizes it so that the DOI is always a link (starts with https://doi.org)
    3. Returns the newly normalized DOI
    """
    doi_url = doi_node["link_url"]
    parsed = urlparse(doi_url)
    if not parsed.scheme and not parsed.netloc:
        parsed_doi = parsed.path
        if not DOI_REGEX_COMPILED.match(parsed_doi):
            errors.append({"link_type": ProjectLinkType.DOI.name, "reason": "Invalid DOI"})
            return None
        doi_url = f"https://doi.org/{parsed_doi}"
    return doi_url


def post(body: dict, user: str):
    errors = []
    doi_url = None
    if doi_node := get_doi_link_node(body, errors):
        if doi_url := portal_get_normalized_doi_url(doi_node, errors):
            doi_node["link_url"] = doi_url
    collection_id = create_collection_common(body, user, doi_url, errors)
    return make_response(jsonify({"collection_id": collection_id}), 201)


@dbconnect
def create_collection_common(body: dict, user: str, doi: str, errors: list):
    db_session = g.db_session
    verify_collection_body(body, errors)
    if doi is not None:
        publisher_metadata = get_publisher_metadata(doi, errors)
    else:
        publisher_metadata = None

    if errors:
        raise InvalidParametersHTTPException(detail=errors)

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

    return collection.id


def get_collection_and_verify_body(db_session: Session, collection_id: str, body: dict, token_info: dict):
    errors = []
    verify_collection_body(body, errors)
    collection = get_collection_else_forbidden(
        db_session,
        collection_id,
        visibility=CollectionVisibility.PRIVATE.name,
        owner=owner_or_allowed(token_info),
    )
    return collection, errors
