from flask import make_response, jsonify

from ....common.entities import Project


def get_projects_list(query_user_uuid: str = '', from_date: int = None, to_date: int = None):
    projects = [p.to_dict() for p in Project.list_in_time_range(from_date, to_date)]
    return make_response(jsonify(projects), 200)


def get_project_details(path_project_uuid: str):
    raise NotImplementedError


def delete_project(path_project_uuid: str):
    raise NotImplementedError


def get_project_dataset(path_dataset_uuid: str):
    raise NotImplementedError
