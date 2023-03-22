from flask import jsonify

from backend.common.utils.http_exceptions import ForbiddenHTTPException, NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import (
    reshape_for_curation_api,
    validate_uuid_else_forbidden,
)
from backend.layers.common.entities import CollectionVersionId
from backend.portal.api.providers import get_business_logic


def get(collection_version_id: str):
    """
    Fetch the full metadata of a published CollectionVersion, including its datasets
    """
    validate_uuid_else_forbidden(collection_version_id)
    version = get_business_logic().get_collection_version(CollectionVersionId(collection_version_id))
    if version is None or version.canonical_collection.tombstoned is True:
        raise NotFoundHTTPException()
    if version.published_at is None:
        raise ForbiddenHTTPException()
    response = reshape_for_curation_api(version, reshape_for_version_endpoint=True)
    return jsonify(response)
