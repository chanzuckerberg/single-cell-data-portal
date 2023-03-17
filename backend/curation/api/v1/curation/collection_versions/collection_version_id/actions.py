from uuid import UUID

from flask import jsonify, make_response

from backend.common.utils.http_exceptions import NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import reshape_for_curation_api
from backend.layers.auth.user_info import UserInfo
from backend.portal.api.providers import get_business_logic


def get(collection_version_id: str, token_info: dict):
    """
    Fetch the full metadata of a CollectionVersion, including its datasets
    """
    collection_version = get_collection_version_else_forbidden(collection_version_id)
    user_info = UserInfo(token_info)
    response = reshape_for_curation_api(collection_version, user_info, use_dataset_version_explorer_urls=True)
    return jsonify(response)
