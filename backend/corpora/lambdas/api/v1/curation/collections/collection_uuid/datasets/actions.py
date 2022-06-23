import re
from flask import g, make_response

from backend.corpora.api_server.db import dbconnect
from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.http_exceptions import InvalidParametersHTTPException, ConflictException
from backend.corpora.common.utils.regex import DATASET_ID_REGEX, CURATOR_TAG_PREFIX_REGEX, EXTENSION_REGEX
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed
from backend.corpora.lambdas.api.v1.common import get_dataset_else_error, get_collection_else_forbidden

REGEX = f"^({DATASET_ID_REGEX}|{CURATOR_TAG_PREFIX_REGEX})\\.{EXTENSION_REGEX}$"


def validate_curator_tag(curator_tag):
    matched = re.match(REGEX, curator_tag)
    return True if matched and matched.groupdict().get("tag") else False


@dbconnect
def patch(token_info: dict, collection_uuid: str, body: dict, curator_tag: str = None, dataset_uuid: str = None):
    db_session = g.db_session
    dataset = get_dataset_else_error(db_session, dataset_uuid, collection_uuid, curator_tag)
    get_collection_else_forbidden(
        db_session, collection_uuid, visibility=CollectionVisibility.PRIVATE, owner=owner_or_allowed(token_info)
    )
    tag = body["curator_tag"]
    if not validate_curator_tag(tag):
        raise InvalidParametersHTTPException("Invalid Curator Tag")

    # Check if the curator_tag is unique across datasets in the collection.
    if dataset.curator_tag != tag:
        if conflict := Dataset.get(db_session, collection_uuid=collection_uuid, curator_tag=tag):
            raise ConflictException(
                f"Curator_tags must be unique within a collection. Dataset={conflict.id} is using curator_tag={tag}"
            )
        else:
            dataset.update(curator_tag=tag)
    return make_response("", 204)
