from typing import List

from flask import jsonify, make_response

from backend.common.utils.http_exceptions import NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import reshape_for_curation_api, validate_uuid_else_forbidden
from backend.layers.common.entities import CollectionId, CollectionVersionWithPublishedDatasets
from backend.layers.common.helpers import get_dataset_versions_with_published_at_and_collection_version_id
from backend.portal.api.providers import get_business_logic


def get(collection_id: str):
    """
    Return all published versions for a Collection
    """
    validate_uuid_else_forbidden(collection_id)

    all_collection_versions = list(
        get_business_logic().get_all_published_collection_versions_from_canonical(CollectionId(collection_id))
    )
    # Determine published_at and collection_version_id for constituent DatasetVersions
    all_collection_versions_with_calculated_values: List[CollectionVersionWithPublishedDatasets] = []
    for collection_version in all_collection_versions:
        collection_version.datasets = get_dataset_versions_with_published_at_and_collection_version_id(
            collection_version.datasets, all_collection_versions
        )  # hack to allow unpacking via **vars() below
        all_collection_versions_with_calculated_values.append(
            CollectionVersionWithPublishedDatasets(**vars(collection_version))
        )

    collection_versions = sorted(
        [
            reshape_for_curation_api(c_v, reshape_for_version_endpoint=True)
            for c_v in all_collection_versions_with_calculated_values
        ],
        key=lambda c: c["published_at"],
        reverse=True,
    )

    if not collection_versions:
        raise NotFoundHTTPException("Collection not found!")

    return make_response(jsonify(collection_versions))
