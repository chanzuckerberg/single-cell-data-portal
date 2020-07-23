from flask import make_response, jsonify

from ....common.utils.exceptions import AuthorizationError
from ....common.entities import Project


def get_submissions_list(query_user_uuid: str = ""):
    results = Project.list_submissions()
    for result in results:
        result["owner_id"] = result.pop("owner")
    result = dict(submissions=results)
    return make_response(jsonify(result), 200)


def create_new_submission(request_body: dict):
    raise NotImplementedError


def get_submission_details(project_uuid: str):
    result = Project.get_project(project_uuid)
    if result:
        result = result.to_dict()
        # Reshape the data to match.
        result["s3_bucket_key"] = result.pop("s3_bucket", None)
        result["owner"] = result.pop("user")
        result["links"] = [dict(url=link["link_url"], type=link["link_type"]) for link in result["links"]]
        result["attestation"] = dict(needed=result.pop("needs_attestation", None), tc_uri=result.pop("tc_uri", None))
        for dataset in result["datasets"]:
            dataset["dataset_deployments"] = dataset.pop("deployment_directories")
            dataset["dataset_assets"] = dataset.pop("artifacts")
            dataset["preprint_doi"] = dict(title=dataset.pop("preprint_doi"))
            dataset["publication_doi"] = dict(title=dataset.pop("publication_doi"))

        return make_response(jsonify(result), 200)
    else:
        raise AuthorizationError()


def delete_submission(path_project_uuid: str):
    raise NotImplementedError


def add_file_to_submission(path_project_uuid: str, request_body: dict):
    raise NotImplementedError


def delete_dataset_from_submission(path_project_uuid: str, path_dataset_uuid: str):
    raise NotImplementedError


def validate_submission(path_project_uuid: str):
    raise NotImplementedError


def save_submission(path_project_uuid: str, request_body: dict):
    raise NotImplementedError


def publish_submission(path_project_uuid: str):
    raise NotImplementedError
