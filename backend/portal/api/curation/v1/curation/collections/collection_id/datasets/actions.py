from flask import make_response, jsonify

from backend.common.utils.http_exceptions import MethodNotAllowedException
from backend.layers.api.providers import get_business_logic
from backend.layers.auth.user_info import UserInfo
from backend.portal.api.curation.v1.curation.collections.common import (
    get_infered_collection_version_else_forbidden,
    is_owner_or_allowed_else_forbidden,
)
from backend.layers.common.entities import CollectionVersionId


def post(token_info: dict, collection_id: str):
    user_info = UserInfo(token_info)
    business_logic = get_business_logic()

    collection_version = get_infered_collection_version_else_forbidden(collection_id)
    is_owner_or_allowed_else_forbidden(collection_version, user_info)

    if collection_version.published_at is not None:
        raise MethodNotAllowedException("Collection must be PRIVATE Collection, or a revision of a PUBLIC Collection.")

    dataset_id, _ = business_logic.create_empty_dataset(CollectionVersionId(collection_id))

    return make_response(jsonify({"id": dataset_id.id}), 201)
