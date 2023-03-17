from flask import jsonify

from backend.curation.api.v1.curation.collections.common import (
    get_collection_version_else_forbidden,
    reshape_for_curation_api,
)


def get(collection_version_id: str):
    """
    Fetch the full metadata of a CollectionVersion, including its datasets
    """
    collection_version = get_collection_version_else_forbidden(collection_version_id)
    response = reshape_for_curation_api(collection_version, use_dataset_version_explorer_urls=True)
    return jsonify(response)
