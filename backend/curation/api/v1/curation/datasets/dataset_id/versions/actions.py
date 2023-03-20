from uuid import UUID

from flask import jsonify, make_response

from backend.common.utils.http_exceptions import NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import reshape_datasets_for_curation_api
from backend.layers.common.entities import DatasetId
from backend.portal.api.providers import get_business_logic


def get(dataset_id: str):
    """
    List all versions for a given Dataset that have been published.
    """
    business_logic = get_business_logic()
    try:
        UUID(dataset_id)
    except ValueError as e:
        raise NotFoundHTTPException("Dataset not found!") from e

    dataset_versions = business_logic.get_prior_published_versions_for_dataset(DatasetId(dataset_id))
    if len(dataset_versions) == 0:
        raise NotFoundHTTPException("Dataset not found!")

    # business function returns in chronological order; must return in reverse chronological
    dataset_versions.reverse()
    response_body = {
        "datasets": reshape_datasets_for_curation_api(dataset_versions, use_canonical_url=False),
    }
    return make_response(jsonify(response_body), 200)
