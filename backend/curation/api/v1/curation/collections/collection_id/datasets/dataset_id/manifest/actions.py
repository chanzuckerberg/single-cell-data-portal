from flask import Response

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
)
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.exceptions import (
    CollectionIsPublishedException,
    CollectionNotFoundException,
    DatasetInWrongStatusException,
    DatasetNotFoundException,
    InvalidURIException,
)
from backend.portal.api.providers import get_business_logic


def put(collection_id: str, dataset_id: str, body: dict, token_info: dict):
    # TODO: deduplicate from ApiCommon. We need to settle the class/module level debate before can do that
    business_logic = get_business_logic()

    collection_version, dataset_version = _get_collection_and_dataset(collection_id, dataset_id)

    if not UserInfo(token_info).is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException()

    try:
        business_logic.ingest_dataset(
            collection_version.version_id,
            body,
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
