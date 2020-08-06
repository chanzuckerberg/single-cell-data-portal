from flask import make_response, jsonify

from ....common.utils.db_utils import db_session
from ....common.utils.exceptions import ForbiddenHTTPException
from ....common.entities import Project


@db_session
def get_submissions_list(query_user_uuid: str = ""):
    results = Project.list_submissions()
    for result in results:
        result["owner_id"] = result.pop("owner")
    result = dict(submissions=results)
    return make_response(jsonify(result), 200)


@db_session
def create_new_submission(request_body: dict):
    raise NotImplementedError


@db_session
def get_submission_details(project_uuid: str):
    project = Project.get_submission(project_uuid)
    if project:
        result = project.reshape_for_api()
        return make_response(jsonify(result), 200)
    else:
        raise ForbiddenHTTPException()


@db_session
def delete_submission(project_uuid: str):
    raise NotImplementedError


@db_session
def add_file_to_submission(project_uuid: str, request_body: dict):
    raise NotImplementedError


@db_session
def delete_dataset_from_submission(project_uuid: str, dataset_uuid: str):
    raise NotImplementedError


@db_session
def validate_submission(project_uuid: str):
    raise NotImplementedError


@db_session
def save_submission(project_uuid: str, request_body: dict):
    raise NotImplementedError


@db_session
def publish_submission(project_uuid: str):
    raise NotImplementedError
