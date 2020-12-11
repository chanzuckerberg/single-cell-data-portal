from flask import make_response, jsonify

from ....common.utils.db_utils import db_session
from ....common.utils.exceptions import ForbiddenHTTPException
from ....common.entities import Collection


@db_session
def get_submissions_list(query_user_uuid: str = ""):
    results = Collection.list_submissions()
    for result in results:
        result["owner_id"] = result.pop("owner")
    result = dict(submissions=results)
    return make_response(jsonify(result), 200)


@db_session
def create_new_submission(request_body: dict):
    raise NotImplementedError


@db_session
def get_submission_details(collection_uuid: str):
    collection = Collection.get_submission(collection_uuid)
    if collection:
        result = collection.reshape_for_api()
        return make_response(jsonify(result), 200)
    else:
        raise ForbiddenHTTPException()


@db_session
def delete_submission(collection_uuid: str):
    raise NotImplementedError


@db_session
def add_file_to_submission(collection_uuid: str, request_body: dict):
    raise NotImplementedError


@db_session
def delete_dataset_from_submission(collection_uuid: str, dataset_uuid: str):
    raise NotImplementedError


@db_session
def validate_submission(collection_uuid: str):
    raise NotImplementedError


@db_session
def save_submission(collection_uuid: str, request_body: dict):
    raise NotImplementedError


@db_session
def publish_submission(collection_uuid: str):
    raise NotImplementedError
