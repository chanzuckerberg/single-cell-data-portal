from typing import Tuple

from flask import Response, jsonify, make_response

from backend.common.utils.exceptions import MaxFileSizeExceededException
from backend.common.utils.http_exceptions import (
    ForbiddenHTTPException,
    GoneHTTPException,
    InvalidParametersHTTPException,
    MethodNotAllowedException,
    NotFoundHTTPException,
    TooLargeHTTPException,
)
from backend.curation.api.v1.curation.collections.common import (
    get_inferred_collection_version,
    reshape_dataset_for_curation_api,
    validate_uuid_else_forbidden,
)
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.exceptions import (
    CollectionIsPublishedException,
    CollectionNotFoundException,
    CollectionUpdateException,
    DatasetInWrongStatusException,
    DatasetIsTombstonedException,
    DatasetNotFoundException,
    InvalidURIException,
)
from backend.layers.common.entities import (
    CollectionVersionWithDatasets,
    DatasetId,
    DatasetVersion,
)
from backend.portal.api.providers import get_business_logic


def get(collection_id: str, dataset_id: str = None):
    collection_version, dataset_version = _get_collection_and_dataset(collection_id, dataset_id)

    # A canonical url should be only used in two cases:
    # 1. the collection version is unpublished but it's not a revision
    # 2. the collection version is published
    use_canonical_url = collection_version.is_initial_unpublished_version() or collection_version.is_published()

    response_body = reshape_dataset_for_curation_api(dataset_version, use_canonical_url)
    return make_response(jsonify(response_body), 200)


def _get_collection_and_dataset(
    collection_id: str, dataset_id: str
) -> Tuple[CollectionVersionWithDatasets, DatasetVersion]:
    """
    Get collection and dataset by their ids. Will look up collection by version and canonical id, and dataset by
    canonical only
    """
    validate_uuid_else_forbidden(collection_id)
    validate_uuid_else_forbidden(dataset_id)
    collection_version = get_inferred_collection_version(collection_id)

    # Extract the dataset from the dataset list.
    dataset_version = None
    for dataset in collection_version.datasets:
        if dataset.dataset_id.id == dataset_id:
            dataset_version = dataset
            break
        if dataset.version_id.id == dataset_id:
            raise ForbiddenHTTPException from None
    if dataset_version is None:
        try:
            get_business_logic().get_dataset_version_from_canonical(DatasetId(dataset_id), get_tombstoned=True)
        except DatasetIsTombstonedException:
            raise GoneHTTPException() from None
        raise NotFoundHTTPException() from None

    return collection_version, dataset_version


def delete(token_info: dict, collection_id: str, dataset_id: str, delete_published: str):
    business_logic = get_business_logic()
    user_info = UserInfo(token_info)

    collection_version, dataset_version = _get_collection_and_dataset(collection_id, dataset_id)

    if not user_info.is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException("Unauthorized")

    if dataset_version.version_id not in [v.version_id for v in collection_version.datasets]:
        raise ForbiddenHTTPException(f"Dataset {dataset_id} does not belong to a collection")

    delete_published_flag = False
    if user_info.is_cxg_admin() and delete_published == "true" and not collection_version.published_at:
        delete_published_flag = True

    try:
        business_logic.remove_dataset_version(
            collection_version.version_id, dataset_version.version_id, delete_published_flag
        )
    except CollectionUpdateException:
        raise MethodNotAllowedException(detail="Cannot delete a public Dataset") from None
    return Response(status=202)
    # End of duplicate block


def put(collection_id: str, dataset_id: str, body: dict, token_info: dict):
    # TODO: deduplicate from ApiCommon. We need to settle the class/module level debate before can do that
    url = body.get("url", body.get("link"))
    business_logic = get_business_logic()

    collection_version, dataset_version = _get_collection_and_dataset(collection_id, dataset_id)

    if not UserInfo(token_info).is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException()

    try:
        business_logic.ingest_dataset(
            collection_version.version_id,
            url,
            None,
            None if dataset_id is None else dataset_version.version_id,
        )
        return Response(status=202)
    except CollectionNotFoundException:
        raise ForbiddenHTTPException() from None
    except CollectionIsPublishedException:
        raise ForbiddenHTTPException() from None
    except DatasetNotFoundException:
        raise NotFoundHTTPException() from None
    except InvalidURIException:
        raise InvalidParametersHTTPException(detail="The dropbox shared link is invalid.") from None
    except MaxFileSizeExceededException:
        raise TooLargeHTTPException() from None
    except DatasetInWrongStatusException:
        raise MethodNotAllowedException(
            detail="Submission failed. A dataset cannot be updated while a previous update for the same dataset "
            "is in progress. Please cancel the current submission by deleting the dataset, or wait until "
            "the submission has finished processing."
        ) from None
    # End of duplicate block
