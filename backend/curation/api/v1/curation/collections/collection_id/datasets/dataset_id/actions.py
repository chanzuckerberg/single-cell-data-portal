from flask import Response, jsonify, make_response

from backend.common.utils.exceptions import MaxFileSizeExceededException
from backend.common.utils.http_exceptions import (
    ForbiddenHTTPException,
    InvalidParametersHTTPException,
    MethodNotAllowedException,
    NotFoundHTTPException,
    TooLargeHTTPException,
)
from backend.curation.api.v1.curation.collections.common import (
    _get_collection_and_dataset,
    reshape_dataset_for_curation_api,
)
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.exceptions import (
    CollectionIsPublishedException,
    CollectionNotFoundException,
    CollectionUpdateException,
    DatasetInWrongStatusException,
    DatasetIsPrivateException,
    DatasetNotFoundException,
    InvalidIngestionManifestException,
    InvalidMetadataException,
    InvalidURIException,
)
from backend.layers.common.entities import (
    DatasetArtifactMetadataUpdate,
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


def delete(token_info: dict, collection_id: str, dataset_id: str, delete_published: bool = False):
    business_logic = get_business_logic()
    user_info = UserInfo(token_info)

    collection_version, dataset_version = _get_collection_and_dataset(collection_id, dataset_id)

    if not user_info.is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException("Unauthorized")

    if dataset_version.version_id not in [v.version_id for v in collection_version.datasets]:
        raise ForbiddenHTTPException(f"Dataset {dataset_id} does not belong to a collection")

    if user_info.is_cxg_admin() and delete_published and not collection_version.published_at:
        try:
            business_logic.remove_dataset_version(
                collection_version.version_id, dataset_version.version_id, delete_published=True
            )
        except DatasetIsPrivateException:
            raise InvalidParametersHTTPException(
                detail='Query param "delete_published=true" is set but Dataset is not published.'
            ) from None
    else:
        try:
            business_logic.remove_dataset_version(collection_version.version_id, dataset_version.version_id)
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
    except InvalidIngestionManifestException as e:
        raise InvalidParametersHTTPException(detail=e.message) from None
    except DatasetInWrongStatusException:
        raise MethodNotAllowedException(
            detail="Submission failed. A dataset cannot be updated while a previous update for the same dataset "
            "is in progress. Please cancel the current submission by deleting the dataset, or wait until "
            "the submission has finished processing."
        ) from None
    # End of duplicate block


def patch(collection_id: str, dataset_id: str, body: dict, token_info: dict):
    """
    Update a dataset's metadata.
    """

    # Find collection and dataset.
    collection_version, dataset_version = _get_collection_and_dataset(collection_id, dataset_id)

    # Confirm user has permission to update dataset.
    if not UserInfo(token_info).is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException()

    # Create payload and attempt update.
    payload = DatasetArtifactMetadataUpdate(body.get("title"))
    try:
        get_business_logic().update_dataset_artifact_metadata(
            collection_version.version_id, dataset_version.version_id, payload
        )
    except InvalidMetadataException as ex:
        raise InvalidParametersHTTPException(ext=dict(invalid_parameters=ex.errors)) from None
    except CollectionNotFoundException:
        raise NotFoundHTTPException() from None
    except CollectionIsPublishedException:
        raise ForbiddenHTTPException() from None
    except DatasetInWrongStatusException:
        raise MethodNotAllowedException(
            detail="Dataset cannot be updated if processing status is not SUCCESS."
        ) from None

    return Response(status=202)
