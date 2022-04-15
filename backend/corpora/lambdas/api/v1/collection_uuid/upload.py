import requests
from flask import make_response, g

from .....common.corpora_config import CorporaConfig
from .....common.corpora_orm import CollectionVisibility, ProcessingStatus
from .....common import upload_sfn
from .....common.entities import Collection, Dataset
from .....api_server.db import dbconnect
from .....common.utils.dl_sources.url import MissingHeaderException, from_url
from .....common.utils.exceptions import (
    ForbiddenHTTPException,
    InvalidParametersHTTPException,
    TooLargeHTTPException,
    MethodNotAllowedException,
    NotFoundHTTPException,
)
from .....common.utils.math_utils import GB
from ..common import owner_or_allowed


def link(collection_uuid: str, body: dict, user: str):
    dataset_id = upload_from_link(collection_uuid, user, body["url"])
    return make_response({"dataset_uuid": dataset_id}, 202)


def relink(collection_uuid: str, body: dict, user: str):
    dataset_id = upload_from_link(collection_uuid, user, body["url"], body["id"])
    return make_response({"dataset_uuid": dataset_id}, 202)


@dbconnect
def upload_from_link(collection_uuid: str, token_info: dict, url: str, dataset_id: str = None):
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

    if resp.get("size") is not None and resp["size"] > CorporaConfig().upload_max_file_size_gb * GB:
        raise TooLargeHTTPException()
    if resp["name"].rsplit(".")[-1].lower() not in CorporaConfig().upload_file_formats:
        raise InvalidParametersHTTPException("The file referred to by the link is not a support file format.")

    # Get the Collection
    collection = Collection.get_collection(
        db_session,
        collection_uuid,
        visibility=CollectionVisibility.PRIVATE,  # Do not allow changes to public Collections
        owner=owner_or_allowed(token_info),
    )
    if not collection:
        raise ForbiddenHTTPException

    if dataset_id:
        # Update dataset
        dataset = Dataset.get(db_session, dataset_id)
        if collection_uuid == dataset.collection_id:
            if dataset.processing_status.processing_status in [ProcessingStatus.SUCCESS, ProcessingStatus.FAILURE]:
                dataset.reprocess()
            else:
                raise MethodNotAllowedException
        else:
            raise NotFoundHTTPException

    else:
        # Add new dataset
        dataset = Dataset.create(db_session, collection=collection)
    dataset.update(processing_status=dataset.new_processing_status())
    # Start processing link
    upload_sfn.start_upload_sfn(collection_uuid, dataset.id, valid_link.url)
    return dataset.id
