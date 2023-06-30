from flask import jsonify, make_response

from backend.common.utils.http_exceptions import GoneHTTPException, NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import (
    reshape_dataset_for_curation_api,
    validate_uuid_else_forbidden,
)
from backend.layers.business.exceptions import DatasetIsTombstonedException
from backend.layers.common.entities import DatasetVersionId
from backend.portal.api.providers import get_business_logic


def get(dataset_version_id: str):
    """
    Fetch Dataset metadata for a Dataset version.
    """
    business_logic = get_business_logic()
    validate_uuid_else_forbidden(dataset_version_id)
    try:
        dataset_version = business_logic.get_prior_published_dataset_version(DatasetVersionId(dataset_version_id))
        if dataset_version is None:
            raise NotFoundHTTPException("Dataset version not found")
    except DatasetIsTombstonedException:
        raise GoneHTTPException() from None

    response_body = reshape_dataset_for_curation_api(dataset_version, use_canonical_url=False, as_canonical=False)
    return make_response(jsonify(response_body), 200)
