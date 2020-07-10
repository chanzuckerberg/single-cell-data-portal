from flask import make_response, jsonify

from ....common.entities import Project


def get_projects_list(user_uuid: str = "", from_date: int = None, to_date: int = None):
    result = dict(projects=Project.list_in_time_range(from_date=from_date, to_date=to_date),)
    if from_date:
        result["from_date"] = from_date
    if to_date:
        result["to_date"] = to_date
    return make_response(jsonify(result), 200)


def get_project_details(project_uuid: str):
    result = Project.get_project(project_uuid).to_dict()

    # Reshape the data to match.
    result["s3_bucket_key"] = result.pop("s3_bucket")
    result["owner"] = result.pop("user")
    result["links"] = [dict(url=link["link_url"], type=link["link_type"]) for link in result["links"]]
    result["attestation"] = dict(needed=result.pop("needs_attestation"), tc_uri=result.pop("tc_uri"))
    for dataset in result["datasets"]:
        dataset["dataset_deployments"] = dataset.pop("deployment_directories")
        dataset["dataset_assets"] = dataset.pop("artifacts")
        dataset["preprint_doi"] = dict(title=dataset.pop("preprint_doi"))
        dataset["publication_doi"] = dict(title=dataset.pop("publication_doi"))

    return make_response(jsonify(result), 200)


def delete_project(path_project_uuid: str):
    raise NotImplementedError


def get_project_dataset(path_dataset_uuid: str):
    raise NotImplementedError
