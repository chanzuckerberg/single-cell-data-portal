from flask import make_response, jsonify

from ....common.utils.exceptions import ForbiddenHTTPException
from ....common.entities import Project, Dataset


def get_submissions_list(query_user_uuid: str = ""):
    results = Project.list_submissions()
    for result in results:
        result["owner_id"] = result.pop("owner")
    result = dict(submissions=results)
    return make_response(jsonify(result), 200)


def create_new_submission(request_body: dict):
    raise NotImplementedError


def get_submission_details(project_uuid: str):
    project = Project.get_submission(project_uuid)
    if project:
        result = project.reshape_for_api()
        return make_response(jsonify(result), 200)
    else:
        raise ForbiddenHTTPException()


def delete_submission(project_uuid: str):
    project = Project.get_submission(project_uuid)
    #  TODO delete uploaded files if they are not associated with a published project
    #   see https://github.com/chanzuckerberg/corpora-data-portal/issues/506
    if project:
        project.delete()
        return make_response("", 202)
    else:
        raise ForbiddenHTTPException()


def add_file_to_submission(project_uuid: str, request_body: dict):
    raise NotImplementedError


def delete_dataset_from_submission(project_uuid: str, dataset_uuid: str):
    dataset = Dataset.get(dataset_uuid)
    # TODO delete uploaded files if they are not assoicated with a published project
    #   see https://github.com/chanzuckerberg/corpora-data-portal/issues/506
    if dataset and dataset.project_id == project_uuid and dataset.is_submission():
        dataset.delete()
        return make_response("", 202)
    else:
        raise ForbiddenHTTPException()


def validate_submission(project_uuid: str):
    raise NotImplementedError


def save_submission(project_uuid: str, request_body: dict):
    raise NotImplementedError


def publish_submission(project_uuid: str):
    raise NotImplementedError
