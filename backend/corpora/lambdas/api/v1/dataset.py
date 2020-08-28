from flask import make_response, jsonify

from ....common.entities import Dataset
from ....common.utils.db_utils import db_session
from ....common.utils.exceptions import ForbiddenHTTPException, NotFoundHTTPException
from ....common.utils.s3_utils import generate_file_url, head_file


@db_session
def post_dataset_artifact(dataset_uuid: str, artifact_uuid):
    dataset = Dataset.get(dataset_uuid)
    if not dataset:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}' not found.")
    artifact = [artifact for artifact in dataset.artifacts if artifact.id == artifact_uuid]
    if not artifact:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}/artifact/{artifact_uuid}' not found.")
    else:
        artifact = artifact[0]
    presigned_url = generate_file_url(artifact.s3_uri)
    return make_response(jsonify(download_link=presigned_url), 201)

@db_session
def head_dataset_artifact(dataset_uuid: str, artifact_uuid):
    dataset = Dataset.get(dataset_uuid)
    if not dataset:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}' not found.")
    artifact = [artifact for artifact in dataset.artifacts if artifact.id == artifact_uuid]
    if not artifact:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}/artifact/{artifact_uuid}' not found.")
    else:
        artifact = artifact[0]
    head_response = head_file(artifact.s3_uri)
    return make_response(jsonify(download_link=head_response), 201)

