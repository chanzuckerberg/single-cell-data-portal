from flask import jsonify, make_response

from backend.curation.api.v1.curation.collections.common import reshape_for_curation_api
from backend.layers.common.entities import CollectionId
from backend.portal.api.providers import get_business_logic


def get(collection_id: str):
    """
    Return all published versions for a Collection
    """
    collection_versions = list(
        map(
            reshape_for_curation_api,
            get_business_logic().get_all_published_collection_versions_from_canonical(CollectionId(collection_id)),
        )
    )

    collection_versions.sort(key=lambda c: c["published_at"])

    return make_response(jsonify(collection_versions))
