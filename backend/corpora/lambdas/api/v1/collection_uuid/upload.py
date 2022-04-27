import requests
from flask import make_response, g

from .....api_server.db import dbconnect
from .....common.upload_sfn import upload
from .....common.utils.dl_sources.url import MissingHeaderException, from_url
from .....common.utils.exceptions import (
    ForbiddenHTTPException,
    InvalidParametersHTTPException,
    TooLargeHTTPException,
    MethodNotAllowedException,
    NotFoundHTTPException,
    MaxFileSizeExceededException,
    InvalidFileFormatException,
    NonExistentCollectionException,
    InvalidProcessingStateException,
    NonExistentDatasetException,
)


def link(collection_uuid: str, body: dict, user: str):
    dataset_id = upload_from_link(collection_uuid, user, body["url"])
    return make_response({"dataset_uuid": dataset_id}, 202)


def relink(collection_uuid: str, body: dict, user: str):
    dataset_id = upload_from_link(collection_uuid, user, body["url"], body["id"])
    return make_response({"dataset_uuid": dataset_id}, 202)


@dbconnect
def upload_from_link(collection_uuid: str, user: str, url: str, dataset_id: str = None):
    db_session = g.db_session
    # Verify Dropbox URL
    valid_link = from_url(url)
    if not valid_link:
        raise InvalidParametersHTTPException("The dropbox shared link is invalid.")

    # Get file info
    try:
        resp = valid_link.file_info()
    except requests.HTTPError:
        raise InvalidParametersHTTPException("The URL provided causes an error with Dropbox.")
    except MissingHeaderException as ex:
        raise InvalidParametersHTTPException(ex.detail)

    file_size = resp.get("size")
    file_extension = resp["name"].rsplit(".")[-1].lower()

    try:
        return upload(db_session, collection_uuid, user, url, file_size, file_extension, dataset_id)
    except MaxFileSizeExceededException:
        raise TooLargeHTTPException()
    except InvalidFileFormatException:
        raise InvalidParametersHTTPException("The file referred to by the link is not a support file format.")
    except NonExistentCollectionException:
        raise ForbiddenHTTPException()
    except InvalidProcessingStateException:
        raise MethodNotAllowedException()
    except NonExistentDatasetException:
        raise NotFoundHTTPException()
