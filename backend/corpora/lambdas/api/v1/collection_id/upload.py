import requests
from flask import make_response, g

from .....api_server.db import dbconnect
from .....common.upload import upload
from .....common.utils.dl_sources.url import MissingHeaderException, from_url
from .....common.utils.exceptions import (
    MaxFileSizeExceededException,
    InvalidFileFormatException,
    NonExistentCollectionException,
    InvalidProcessingStateException,
    NonExistentDatasetException,
)
from .....common.utils.http_exceptions import (
    ForbiddenHTTPException,
    InvalidParametersHTTPException,
    TooLargeHTTPException,
    MethodNotAllowedException,
    NotFoundHTTPException,
)


def link(collection_id: str, body: dict, token_info: dict):
    dataset_id = upload_from_link(collection_id, token_info, body["url"], curator_tag=body.get("curator_tag"))
    return make_response({"dataset_id": dataset_id}, 202)


def relink(collection_id: str, body: dict, token_info: dict):
    dataset_id = upload_from_link(
        collection_id,
        token_info,
        body.get("url", body.get("link")),
        body.get("id"),
        curator_tag=body.get("curator_tag"),
    )
    return make_response({"dataset_id": dataset_id}, 202)


@dbconnect
def upload_from_link(collection_id: str, token_info: dict, url: str, dataset_id: str = None, curator_tag: str = None):
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
    file_extension = resp["name"].rsplit(".")[-1].lower()

    try:
        return upload(
            db_session,
            collection_id=collection_id,
            url=url,
            file_size=file_size,
            file_extension=file_extension,
            user=token_info["sub"],
            scope=token_info["scope"],
            dataset_id=dataset_id,
            curator_tag=curator_tag,
        )
    except MaxFileSizeExceededException:
        raise TooLargeHTTPException()
    except InvalidFileFormatException:
        raise InvalidParametersHTTPException(detail="The file referred to by the link is not a support file format.")
    except NonExistentCollectionException:
        raise ForbiddenHTTPException()
    except InvalidProcessingStateException:
        raise MethodNotAllowedException()
    except NonExistentDatasetException:
        raise NotFoundHTTPException()
