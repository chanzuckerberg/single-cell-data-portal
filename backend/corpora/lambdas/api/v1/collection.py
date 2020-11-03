from flask import make_response, jsonify

from ....common.utils.db_utils import db_session
from ....common.entities import Collection
from ....common.utils.exceptions import ForbiddenHTTPException


@db_session
def get_collections_list(user_uuid: str = "", from_date: int = None, to_date: int = None):
    result = dict(collections=Collection.list_collections_in_time_range(from_date=from_date, to_date=to_date))
    if from_date:
        result["from_date"] = from_date
    if to_date:
        result["to_date"] = to_date
    return make_response(jsonify(result), 200)


@db_session
def get_collection_details(collection_uuid: str):
    collection = Collection.get_collection(collection_uuid)
    if collection:
        result = collection.reshape_for_api()
        return make_response(jsonify(result), 200)
    else:
        raise ForbiddenHTTPException()


@db_session
def delete_collection(collection_uuid: str):
    raise NotImplementedError


@db_session
def get_collection_dataset(dataset_uuid: str):
    raise NotImplementedError
