from flask import jsonify, make_response

from backend.common.utils.http_exceptions import GoneHTTPException, NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import reshape_for_curation_api, validate_uuid_else_forbidden
from backend.layers.common.entities import CollectionId
from backend.portal.api.providers import get_business_logic


def get(collection_id: str):
    """
    Return all published versions for a Collection
    """
    validate_uuid_else_forbidden(collection_id)
    business_logic = get_business_logic()
    canonical_collection = business_logic.get_canonical_collection(CollectionId(collection_id))
    if not canonical_collection:
        raise NotFoundHTTPException("Collection not found!")
    if canonical_collection.tombstoned:
        raise GoneHTTPException()

    collection_versions = sorted(
        [
            reshape_for_curation_api(c_v, reshape_for_version_endpoint=True)
            for c_v in get_business_logic().get_all_published_collection_versions_from_canonical(
                canonical_collection.id
            )
        ],
        key=lambda c: c["published_at"],
        reverse=True,
    )

    if not collection_versions:
        raise NotFoundHTTPException("Collection not found!")

    return make_response(jsonify(collection_versions))
