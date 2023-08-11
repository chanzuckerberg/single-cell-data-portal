from flask import jsonify, make_response

from backend.common.utils.http_exceptions import InvalidParametersHTTPException
from backend.curation.api.v1.curation.collections.common import extract_doi_from_links, reshape_dataset_for_curation_api
from backend.portal.api.providers import get_business_logic


def get(schema_version: str = None):  # type: ignore
    """
    Datasets index endpoint to retrieve full metadata for all public Datasets.
    """
    all_datasets_with_collection_name_and_doi = []
    if not schema_version:
        collections_with_datasets = get_business_logic().get_all_mapped_collection_versions_with_datasets()
    else:
        version_parts = schema_version.split(".")
        if len(version_parts) > 3 or not all(part.isdigit() for part in version_parts):
            raise InvalidParametersHTTPException(detail="Invalid Schema Version Input")
        while len(version_parts) < 3:
            # wildcard match for exactly 1 character
            version_parts.append("_")
        schema_version = ".".join(version_parts)
        collections_with_datasets = get_business_logic().get_latest_published_collection_versions_by_schema(
            schema_version
        )

    for collection in collections_with_datasets:
        doi, _ = extract_doi_from_links(collection.metadata.links)

        collection_info = {
            "collection_id": collection.collection_id.id,
            "collection_name": collection.metadata.name,
            "collection_doi": doi,
        }

        for dataset in collection.datasets:
            dataset_response_obj = reshape_dataset_for_curation_api(dataset, use_canonical_url=True)
            dataset_response_obj.update(collection_info)
            all_datasets_with_collection_name_and_doi.append(dataset_response_obj)

    return make_response(
        jsonify(
            sorted(
                all_datasets_with_collection_name_and_doi,
                key=lambda d: (
                    d["published_at"],
                    d["dataset_id"],
                ),  # Secondary sort by dataset_id for consistency since some Datasets from the same Collection will have identical published_at dates
                reverse=True,
            )
        ),
        200,
    )
