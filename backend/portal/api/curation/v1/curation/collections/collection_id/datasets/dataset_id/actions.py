from typing import Tuple

from flask import Response, jsonify, make_response

from backend.common.utils.exceptions import MaxFileSizeExceededException
from backend.common.utils.http_exceptions import (
    ForbiddenHTTPException,
    InvalidParametersHTTPException,
    MethodNotAllowedException,
    NotFoundHTTPException,
    TooLargeHTTPException,
)
from backend.layers.api.router import get_business_logic
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.business.exceptions import (
    CollectionIsPublishedException,
    CollectionNotFoundException,
    CollectionUpdateException,
    DatasetInWrongStatusException,
    DatasetNotFoundException,
    InvalidURIException,
)
from backend.layers.common.entities import (
    CollectionId,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetVersion,
    DatasetVersionId,
)
from backend.portal.api.curation.v1.curation.collections.common import (
    get_infered_dataset_version,
    reshape_dataset_for_curation_api,
)


def get(collection_id: str, dataset_id: str = None):
    business_logic = get_business_logic()

    collection_version = business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
    if collection_version is None:
        raise NotFoundHTTPException("Collection not found!")

    version = get_infered_dataset_version(dataset_id)
    if version is None:
        raise NotFoundHTTPException("Dataset not found")

    response_body = reshape_dataset_for_curation_api(version)
    return make_response(jsonify(response_body), 200)


def _get_collection_and_dataset(
    business_logic: BusinessLogicInterface, collection_id: CollectionVersionId, dataset_id: DatasetVersionId
) -> Tuple[CollectionVersionWithDatasets, DatasetVersion]:
    dataset_version = business_logic.get_dataset_version(dataset_id)
    if dataset_version is None:
        raise ForbiddenHTTPException()

    collection_version = business_logic.get_collection_version(collection_id)
    if collection_version is None:
        raise ForbiddenHTTPException()

    return collection_version, dataset_version


def delete(token_info: dict, collection_id: str, dataset_id: str = None):
    business_logic = get_business_logic()
    user_info = UserInfo(token_info)

    collection_version, dataset_version = _get_collection_and_dataset(
        business_logic, CollectionVersionId(collection_id), DatasetVersionId(dataset_id)
    )

    if not user_info.is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException("Unauthorized")
    # End of duplicate block

    # TODO: deduplicate from ApiCommon. We need to settle the class/module level debate before can do that
    if dataset_version.version_id not in [v.version_id for v in collection_version.datasets]:
        raise ForbiddenHTTPException(f"Dataset {dataset_id} does not belong to a collection")

    try:
        business_logic.remove_dataset_version(collection_version.version_id, dataset_version.version_id)
    except CollectionUpdateException:
        raise MethodNotAllowedException(detail="Cannot delete a public Dataset")
    return Response(status=202)
    # End of duplicate block


def put(collection_id: str, dataset_id: str, body: dict, token_info: dict):
    # TODO: deduplicate from ApiCommon. We need to settle the class/module level debate before can do that
    url = body.get("url", body.get("link"))
    business_logic = get_business_logic()

    collection_version, _ = _get_collection_and_dataset(
        business_logic, CollectionVersionId(collection_id), DatasetVersionId(dataset_id)
    )

    if not UserInfo(token_info).is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException()

    try:
        business_logic.ingest_dataset(
            collection_version.version_id,
            url,
            None if dataset_id is None else DatasetVersionId(dataset_id),
        )
        return Response(status=202)
    except CollectionNotFoundException:
        raise ForbiddenHTTPException()
    except CollectionIsPublishedException:
        raise ForbiddenHTTPException()
    except DatasetNotFoundException:
        raise NotFoundHTTPException()
    except InvalidURIException:
        raise InvalidParametersHTTPException(detail="The dropbox shared link is invalid.")
    except MaxFileSizeExceededException:
        raise TooLargeHTTPException()
    except DatasetInWrongStatusException:
        raise MethodNotAllowedException(
            detail="Submission failed. A dataset cannot be updated while a previous update for the same dataset "
            "is in progress. Please cancel the current submission by deleting the dataset, or wait until "
            "the submission has finished processing."
        )
    # End of duplicate block
