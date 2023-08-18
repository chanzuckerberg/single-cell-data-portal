from flask import jsonify, make_response

from backend.common.utils.http_exceptions import GoneHTTPException, NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import (
    reshape_datasets_for_curation_api,
    validate_uuid_else_forbidden,
)
from backend.layers.business.exceptions import DatasetIsTombstonedException
from backend.layers.common.entities import DatasetId
from backend.portal.api.providers import get_business_logic


def get(dataset_id: str):
    """
    List all versions for a given Dataset that have been published.
    """
    business_logic = get_business_logic()
    validate_uuid_else_forbidden(dataset_id)
    try:
        dataset_versions = business_logic.get_prior_published_versions_for_dataset(DatasetId(dataset_id))
    except DatasetIsTombstonedException:
        raise GoneHTTPException() from None
    if len(dataset_versions) == 0:
        raise NotFoundHTTPException("Dataset not found!")

    # business function returns in chronological order; must return in reverse chronological
    dataset_versions.reverse()

    return make_response(
        jsonify(reshape_datasets_for_curation_api(dataset_versions, use_canonical_url=False, as_version=True)), 200
    )
