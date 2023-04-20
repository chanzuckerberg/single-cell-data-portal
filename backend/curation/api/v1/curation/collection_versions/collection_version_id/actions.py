from flask import jsonify

from backend.common.utils.http_exceptions import GoneHTTPException, NotFoundHTTPException
from backend.curation.api.v1.curation.collections.common import (
    reshape_for_curation_api,
    validate_uuid_else_forbidden,
)
from backend.layers.business.exceptions import (
    CollectionIsTombstonedException,
    PublishedCollectionVersionNotFoundException,
)
from backend.layers.common.entities import CollectionVersionId
from backend.portal.api.providers import get_business_logic


def get(collection_version_id: str):
    """
    Fetch the full metadata of a published CollectionVersion, including its datasets
    """
    validate_uuid_else_forbidden(collection_version_id)
    try:
        version = get_business_logic().get_published_collection_version__discover_api(
            CollectionVersionId(collection_version_id)
        )
    except PublishedCollectionVersionNotFoundException as e:
        raise NotFoundHTTPException() from e
    except CollectionIsTombstonedException as e:
        raise GoneHTTPException() from e
    response = reshape_for_curation_api(version, reshape_for_version_endpoint=True)
    return jsonify(response)
