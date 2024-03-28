from flask import jsonify, make_response

from backend.common.utils.http_exceptions import ForbiddenHTTPException, InvalidParametersHTTPException
from backend.curation.api.v1.curation.collections.common import reshape_dataset_for_curation_datasets_index_api
from backend.layers.auth.user_info import UserInfo
from backend.layers.common.entities import DatasetVisibility
from backend.portal.api.providers import get_business_logic


def get(token_info: dict, schema_version: str = None, visibility: str = None):
    """
    Datasets index endpoint to retrieve full metadata. Only return Dataset data for which the curator is authorized.
    :param token_info: access token info.
    :param schema_version: the schema version to filter the datasets by, PUBLIC Datasets only.
    :param visibility: the DatasetVisibility in string form.
    """

    # Handle retrieval of private datasets.
    if visibility == DatasetVisibility.PRIVATE.name:
        if schema_version:
            raise InvalidParametersHTTPException(detail="schema_version is not allowed for PRIVATE Datasets.")

        user_info = UserInfo(token_info)
        if user_info.is_none():
            raise ForbiddenHTTPException(detail="Not authorized to query for PRIVATE Dataset.")

        owner = None
        if not user_info.is_super_curator():  # No owner if user is super curator.
            owner = user_info.user_id

        collections_with_datasets = get_business_logic().get_private_collection_versions_with_datasets(owner)
    # Handle retrieval of public datasets.
    else:
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

    # Shape datasets for response.
    all_datasets_with_collection_name_and_doi = []
    for collection in collections_with_datasets:
        for dataset in collection.datasets:
            dataset_response_obj = reshape_dataset_for_curation_datasets_index_api(visibility, collection, dataset)
            all_datasets_with_collection_name_and_doi.append(dataset_response_obj)

    return make_response(
        jsonify(
            sorted(
                all_datasets_with_collection_name_and_doi,
                key=lambda d: (
                    d["published_at"] is None,
                    d["published_at"],
                    d["dataset_id"],
                ),  # Secondary sort by dataset_id for consistency since some Datasets from the same Collection will have identical published_at dates
                reverse=True,
            )
        ),
        200,
    )
