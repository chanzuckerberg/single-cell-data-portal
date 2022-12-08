from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.corpora_orm import CollectionVisibility, ProcessingStatus
from backend.common.entities import Dataset
from backend.common.utils.http_exceptions import MethodNotAllowedException
from backend.layers.api.router import get_business_logic
from backend.layers.auth.user_info import UserInfo
from backend.portal.api.app.v1.authorization import owner_or_allowed
from backend.portal.api.collections_common import get_collection_else_forbidden
from backend.portal.api.curation.v1.curation.collections.common import get_infered_collection_version_else_forbidden, is_owner_or_allowed_else_forbidden


@dbconnect
def post(token_info: dict, collection_id: str):
    user_info = UserInfo(token_info)
    business_logic = get_business_logic()
    collection_version = get_infered_collection_version_else_forbidden(collection_id)
    is_owner_or_allowed_else_forbidden(collection_version, user_info)

    if collection_version.published_at is not None:
        raise MethodNotAllowedException("Collection must be PRIVATE Collection, or a revision of a PUBLIC Collection.")

    business_logic.create_da

    dataset = Dataset.create(
        db_session, collection=collection, processing_status={"processing_status": ProcessingStatus.INITIALIZED}
    )
    return make_response(jsonify({"id": dataset.id}), 201)
