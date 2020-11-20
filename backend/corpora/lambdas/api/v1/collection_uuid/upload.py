from flask import make_response, jsonify
from typing import Optional

from .....common.corpora_orm import DbCollection, CollectionVisibility
from .....common.utils.db_utils import db_session
from .....common.entities import Collection, Dataset
from .....common.utils.exceptions import ForbiddenHTTPException

@db_session
def link(collection_uuid: str, body: dict, user: str):
    collection = Collection.if_owner(collection_uuid, CollectionVisibility.PRIVATE, user)
    if collection:
        url = body['url']
        dataset = Dataset.create(processing_status=Dataset.new_processing_status(),
                                 collection=collection)
        return make_response({'dataset_uuid': dataset.id}, 202)
    else:
        raise ForbiddenHTTPException
