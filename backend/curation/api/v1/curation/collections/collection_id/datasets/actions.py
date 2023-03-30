from flask import jsonify, make_response

from backend.common.utils.http_exceptions import MethodNotAllowedException
from backend.curation.api.v1.curation.collections.common import (
    get_inferred_collection_version,
    is_owner_or_allowed_else_forbidden,
)
from backend.layers.auth.user_info import UserInfo
from backend.portal.api.providers import get_business_logic


def post(token_info: dict, collection_id: str):
    user_info = UserInfo(token_info)
    business_logic = get_business_logic()

    collection_version = get_inferred_collection_version(collection_id)
    is_owner_or_allowed_else_forbidden(collection_version, user_info)

    if collection_version.published_at is not None:
        raise MethodNotAllowedException("Collection must be PRIVATE Collection, or a revision of a PUBLIC Collection.")

    _, dataset_id = business_logic.create_empty_dataset(collection_version.version_id)

    return make_response(jsonify({"dataset_id": dataset_id.id}), 201)
