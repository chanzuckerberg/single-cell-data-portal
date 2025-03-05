from flask import jsonify, make_response

from backend.common.utils.http_exceptions import GoneHTTPException, NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import (
    validate_uuid_else_forbidden,
)
from backend.layers.business.exceptions import (
    ArtifactNotFoundException,
    DatasetIsTombstonedException,
    DatasetNotFoundException,
)
from backend.layers.common.entities import DatasetVersionId
from backend.portal.api.providers import get_business_logic


def get(dataset_version_id: str):
    """
    Fetch anndata metadata for a Dataset version.
    """
    business_logic = get_business_logic()
    validate_uuid_else_forbidden(dataset_version_id)
    try:
        anndata_metadata = business_logic.get_anndata_metadata(DatasetVersionId(dataset_version_id))
    except DatasetNotFoundException:
        raise NotFoundHTTPException("Dataset version not found") from None
    except ArtifactNotFoundException:
        raise NotFoundHTTPException("H5AD artifact not found") from None
    except DatasetIsTombstonedException:
        raise GoneHTTPException() from None

    response = {
        "obs_column_names": anndata_metadata.obs_column_names,
        "var_column_names": anndata_metadata.var_column_names,
    }
    return make_response(jsonify(response), 200)
