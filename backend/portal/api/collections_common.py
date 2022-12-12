import logging
import re
from typing import Optional
from urllib.parse import urlparse

import sqlalchemy
from flask import g
from sqlalchemy.orm import Session

from backend.api_server.db import dbconnect
from backend.common.corpora_orm import ProjectLinkType, CollectionVisibility
from backend.common.entities import Dataset, Collection
from backend.common.providers.crossref_provider import CrossrefProvider, CrossrefDOINotFoundException, CrossrefException
from backend.common.utils.http_exceptions import (
    ForbiddenHTTPException,
    MethodNotAllowedException,
    InvalidParametersHTTPException,
    NotFoundHTTPException,
)
from backend.common.utils.regex import CONTROL_CHARS
from backend.portal.api.app.v1.authorization import owner_or_allowed


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


def get_collection_else_forbidden(db_session, collection_id, **kwargs):
    collection = Collection.get_collection(db_session, collection_id, **kwargs)
    if not collection:
        raise ForbiddenHTTPException()
    return collection


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


def delete_dataset_common(db_session: Session, dataset: Dataset, token_info: dict):
    if not dataset:
        raise ForbiddenHTTPException()
    get_collection_else_forbidden(
        db_session,
        dataset.collection.id,
        owner=owner_or_allowed(token_info),
    )
    if dataset.collection.visibility == CollectionVisibility.PUBLIC:
        raise MethodNotAllowedException(detail="Cannot delete a public Dataset")
    if dataset.tombstone is False:
        if dataset.published:
            dataset.update(tombstone=True, published=False)
        else:
            if dataset.original_id:
                # The dataset is a revision of a published dataset
                original = Dataset.get(db_session, dataset.original_id)
                original.create_revision(dataset.collection.id)  # Restore the original dataset and S3 assets
            dataset.asset_deletion()  # Delete the S3 assets and database rows.
            dataset.delete()  # Delete the dataset row.


def get_dataset_else_error(db_session, dataset_id, collection_id, **kwargs) -> Dataset:
    try:
        dataset = Dataset.get(db_session, dataset_id, **kwargs)
    except ValueError:
        raise InvalidParametersHTTPException()
    if not dataset:
        get_collection_else_forbidden(db_session, collection_id)  # if dataset not found, check if the collection exists
        raise NotFoundHTTPException(detail="Dataset not found.")
    return dataset


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
        consortia=body.get("consortia", []),
    )

    return collection.id


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
email_regex = re.compile(r"(.+)@(.+)\.(.+)")


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
