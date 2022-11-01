import requests
from flask import make_response, g

from backend.api_server.db import dbconnect
from backend.common.upload import upload
from backend.common.utils.dl_sources.url import MissingHeaderException, from_url
from backend.common.utils.exceptions import (
    MaxFileSizeExceededException,
    InvalidFileFormatException,
    NonExistentCollectionException,
    InvalidProcessingStateException,
    NonExistentDatasetException,
)
from backend.common.utils.http_exceptions import (
    ForbiddenHTTPException,
    InvalidParametersHTTPException,
    TooLargeHTTPException,
    MethodNotAllowedException,
    NotFoundHTTPException,
)


@dbconnect
def upload_from_link(collection_id: str, token_info: dict, url: str, dataset_id: str = None):
    db_session = g.db_session

    # Verify Dropbox URL
    valid_link = from_url(url)
    if not valid_link:
        raise InvalidParametersHTTPException(detail="The dropbox shared link is invalid.")

    # Get file info
    try:
        resp = valid_link.file_info()
    except requests.HTTPError:
        raise InvalidParametersHTTPException(detail="The URL provided causes an error with Dropbox.")
    except MissingHeaderException as ex:
        raise InvalidParametersHTTPException(detail=ex.detail)

    file_size = resp.get("size")

    try:
        return upload(
            db_session,
            collection_id=collection_id,
            url=url,
            file_size=file_size,
            user=token_info["sub"],
            scope=token_info["scope"],
            dataset_id=dataset_id,
        )
    except MaxFileSizeExceededException:
        raise TooLargeHTTPException()
    except InvalidFileFormatException:
        raise InvalidParametersHTTPException(detail="The file referred to by the link is not a support file format.")
    except NonExistentCollectionException:
        raise ForbiddenHTTPException()
    except InvalidProcessingStateException:
        raise MethodNotAllowedException(
            detail="Submission failed. A dataset cannot be updated while a previous update for the same dataset is in "
            "progress. Please cancel the current submission by deleting the dataset, or wait until the submission has "
            "finished processing.",
        )
    except NonExistentDatasetException:
        raise NotFoundHTTPException()


def post(collection_id: str, body: dict, token_info: dict):
    dataset_id = upload_from_link(collection_id, token_info, body["url"])
    return make_response({"dataset_id": dataset_id}, 202)


def put(collection_id: str, body: dict, token_info: dict):
    dataset_id = upload_from_link(
        collection_id,
        token_info,
        body.get("url", body.get("link")),
        body.get("id"),
    )
    return make_response({"dataset_id": dataset_id}, 202)
