from flask import jsonify, make_response

from backend.curation.api.v1.curation.collections.common import reshape_dataset_for_curation_api_as_canonical
from backend.layers.common.entities import CollectionLinkType
from backend.portal.api.providers import get_business_logic


def get():
    """
    Datasets index endpoint to retrieve full metadata for all public Datasets.
    """
    collections_with_datasets = get_business_logic().get_all_mapped_collection_versions_with_datasets()

    all_datasets_with_collection_name_and_doi = []

    for collection in collections_with_datasets:

        collection_doi = None
        for link in collection.metadata.links:
            if link.type == CollectionLinkType.DOI.name:
                collection_doi = link.uri

        collection_info = {
            "collection_id": collection.collection_id.id,
            "collection_name": collection.metadata.name,
            "collection_doi": collection_doi,
        }

        for dataset in collection.datasets:
            dataset_response_obj = reshape_dataset_for_curation_api_as_canonical(dataset, use_canonical_url=True)
            dataset_response_obj.update(collection_info)
            all_datasets_with_collection_name_and_doi.append(dataset_response_obj)

    return make_response(
        jsonify(sorted(all_datasets_with_collection_name_and_doi, key=lambda d: d["published_at"], reverse=True)), 200
    )
