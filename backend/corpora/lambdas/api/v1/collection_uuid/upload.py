from flask import make_response

from .....common.corpora_config import CorporaConfig
from .....common.corpora_orm import CollectionVisibility
from .....common import upload_sfn
from .....common.utils.db_utils import db_session
from .....common.entities import Collection, Dataset
from .....common.utils import dropbox
from .....common.utils.exceptions import ForbiddenHTTPException, InvalidParametersHTTPException, TooLargeHTTPException
from .....common.utils.math_utils import GB


@db_session
def link(collection_uuid: str, body: dict, user: str):

    # Verify Dropbox URL
    url = body["url"]
    if not dropbox.verify(url):
        raise InvalidParametersHTTPException("The dropbox shared link is invalid.")

    # Get file info
    url = dropbox.get_download_url_from_shared_link(url)
    resp = dropbox.get_file_info(url)
    if resp["size"] > CorporaConfig().upload_max_file_size * GB:
        raise TooLargeHTTPException()
    if resp["name"].rsplit(".")[-1].lower() not in CorporaConfig().upload_file_formats:
        raise InvalidParametersHTTPException("The file referred to by the link is not a support file format.")

    # Create dataset
    collection = Collection.if_owner(collection_uuid, CollectionVisibility.PRIVATE, user)
    if not collection:
        raise ForbiddenHTTPException
    dataset = Dataset.create(processing_status=Dataset.new_processing_status(), collection=collection)

    # Start processing link
    upload_sfn.start_upload_sfn(collection_uuid, dataset.id, url)
    return make_response({"dataset_uuid": dataset.id}, 202)
