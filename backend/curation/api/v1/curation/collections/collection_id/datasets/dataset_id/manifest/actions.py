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
)
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.business import BusinessLogic
from backend.layers.business.exceptions import (
    CollectionIsPublishedException,
    CollectionNotFoundException,
    DatasetInWrongStatusException,
    DatasetNotFoundException,
    InvalidIngestionManifestException,
    InvalidURIException,
)
from backend.layers.common.entities import (
    DatasetArtifactType,
    DatasetVersion,
)
from backend.portal.api.providers import get_business_logic


def get_single_artifact_permanent_url(
    dataset_version: DatasetVersion, artifact_type: DatasetArtifactType, required: bool = False
) -> str | None:
    """Find exactly one artifact of the given type, then return it's canonical URI.

    If `required` is True and no artifact is found, raises ValueError.
    If more than one is found, always raises ValueError.
    """
    artifacts = dataset_version.artifacts
    matches = [a for a in artifacts if a.type == artifact_type]

    if len(matches) > 1:
        raise ValueError(f"Multiple '{artifact_type}' artifacts found.")

    if not matches and required:
        raise ValueError(f"No '{artifact_type}' artifact found.")

    if matches:
        return BusinessLogic.generate_permanent_url(dataset_version, matches[0].id, artifact_type)


def get(collection_id: str, dataset_id: str = None):
    _, dataset_version = _get_collection_and_dataset(collection_id, dataset_id)

    response_body = {}
    for key, artifact_type in [
        ("anndata", DatasetArtifactType.H5AD),
        ("atac_fragment", DatasetArtifactType.ATAC_FRAGMENT),
    ]:
        if uri := get_single_artifact_permanent_url(dataset_version, artifact_type):
            response_body[key] = uri

    return make_response(jsonify(response_body), 200)


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
    except InvalidURIException as e:
        raise InvalidParametersHTTPException(detail=e.args) from None
    except InvalidIngestionManifestException as e:
        raise InvalidParametersHTTPException(detail=e.message) from None
    except MaxFileSizeExceededException:
        raise TooLargeHTTPException() from None
    except DatasetInWrongStatusException:
        raise MethodNotAllowedException(
            detail="Submission failed. A dataset cannot be updated while a previous update for the same dataset "
            "is in progress. Please cancel the current submission by deleting the dataset, or wait until "
            "the submission has finished processing."
        ) from None
    # End of duplicate block
