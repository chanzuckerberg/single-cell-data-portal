from flask import jsonify, make_response

from backend.curation.api.v1.curation.collections.common import reshape_dataset_for_curation_api
from backend.layers.common.entities import CollectionLinkType
from backend.portal.api.providers import get_business_logic


def get():
    """
    Datasets index endpoint to retrieve full metadata for all public Datasets.
    """
    collections_with_datasets = get_business_logic().get_all_published_collections_with_datasets()

    all_datasets_with_collection_name_and_doi = []

    for collection in collections_with_datasets:
        assert collection.datasets  # Any public Collection should have at least one Dataset
        for dataset in collection.datasets:
            dataset_response_obj = reshape_dataset_for_curation_api(dataset, is_published=True, use_canonical_id=True)
            collection_doi = None
            for link in collection.metadata.links:
                if link.type == CollectionLinkType.DOI.name:
                    collection_doi = link.uri
            dataset_response_obj.update({"collection_name": collection.metadata.name, "collection_doi": collection_doi})
            all_datasets_with_collection_name_and_doi.append(dataset_response_obj)

    return make_response(jsonify(all_datasets_with_collection_name_and_doi), 200)
