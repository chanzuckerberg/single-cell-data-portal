from flask import jsonify

from backend.common.utils.http_exceptions import GoneHTTPException, NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import (
    reshape_for_curation_api,
    validate_uuid_else_forbidden,
)
from backend.layers.common.entities import CollectionVersionId, CollectionVersionWithPublishedDatasets
from backend.layers.common.helpers import get_dataset_versions_with_published_at_and_collection_version_id
from backend.portal.api.providers import get_business_logic


def get(collection_version_id: str):
    """
    Fetch the full metadata of a published CollectionVersion, including its datasets
    """
    validate_uuid_else_forbidden(collection_version_id)
    collection_version = get_business_logic().get_collection_version(
        CollectionVersionId(collection_version_id), get_tombstoned=True
    )
    if collection_version is None or collection_version.published_at is None:
        raise NotFoundHTTPException()
    if collection_version.canonical_collection.tombstoned is True:
        raise GoneHTTPException()

    all_collection_versions = list(
        get_business_logic().get_all_published_collection_versions_from_canonical(
            collection_version.canonical_collection.id
        )
    )
    collection_version.datasets = get_dataset_versions_with_published_at_and_collection_version_id(
        collection_version.datasets, all_collection_versions
    )  # hack to allow unpacking via **vars() below
    response = reshape_for_curation_api(
        CollectionVersionWithPublishedDatasets(**vars(collection_version)), reshape_for_version_endpoint=True
    )
    return jsonify(response)
