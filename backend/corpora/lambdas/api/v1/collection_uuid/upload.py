import requests
from flask import make_response, g

from .....common.corpora_config import CorporaConfig
from .....common.corpora_orm import CollectionVisibility, ProcessingStatus
from .....common import upload_sfn
from .....common.entities import Collection, Dataset
from .....common.utils.dl_sources.url import MissingHeaderException, from_url
from .....common.utils.exceptions import (
    ForbiddenHTTPException,
    InvalidParametersHTTPException,
    TooLargeHTTPException,
    MethodNotAllowedException,
)
from .....common.utils.math_utils import GB


def link(collection_uuid: str, body: dict, user: str):
    db_session = g.db_session
    dataset_id = upload_from_link(db_session, collection_uuid, user, body["url"])
    return make_response({"dataset_uuid": dataset_id}, 202)


def relink(collection_uuid: str, body: dict, user: str):
    db_session = g.db_session
    dataset_id = upload_from_link(db_session, collection_uuid, user, body["url"], body["id"])
    return make_response({"dataset_uuid": dataset_id}, 202)


def upload_from_link(db_session, collection_uuid: str, user: str, url: dict, dataset_id: str = None):
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
    if resp["size"] > CorporaConfig().upload_max_file_size_gb * GB:
        raise TooLargeHTTPException()
    if resp["name"].rsplit(".")[-1].lower() not in CorporaConfig().upload_file_formats:
        raise InvalidParametersHTTPException("The file referred to by the link is not a support file format.")

    # Create dataset
    collection = Collection.get_collection(db_session, collection_uuid, CollectionVisibility.PRIVATE, owner=user)
    if not collection:
        raise ForbiddenHTTPException

    if dataset_id:
        dataset = Dataset.get(dataset_id)
        if dataset.processing_status.processing_status in [ProcessingStatus.SUCCESS, ProcessingStatus.FAILURE]:
            dataset.reprocess()
        else:
            raise MethodNotAllowedException
    else:
        dataset = Dataset.create(db_session, collection=collection)
    dataset.update(processing_status=dataset.new_processing_status())
    # Start processing link
    upload_sfn.start_upload_sfn(collection_uuid, dataset.id, valid_link.url)
    return dataset.id
