from uuid import UUID

from flask import jsonify, make_response

from backend.common.utils.http_exceptions import ForbiddenHTTPException
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.exceptions import CollectionVersionException
from backend.layers.common.entities import CollectionId
from backend.portal.api.providers import get_business_logic


def post(collection_id: str, token_info: dict):
    try:
        UUID(collection_id)
    except ValueError as e:
        raise ForbiddenHTTPException(f"Collection {collection_id} does not exist") from e
    published_collection = get_business_logic().get_published_collection_version(CollectionId(collection_id))

    if published_collection is None:
        raise ForbiddenHTTPException(f"Collection {collection_id} does not exist")

    if not UserInfo(token_info).is_user_owner_or_allowed(published_collection.owner):
        raise ForbiddenHTTPException("Unauthorized")

    try:
        collection_version = get_business_logic().create_collection_version(CollectionId(collection_id))
    except CollectionVersionException:
        raise ForbiddenHTTPException("Another revision is already in progress") from None

    return make_response(jsonify({"collection_id": collection_version.version_id.id}), 201)
