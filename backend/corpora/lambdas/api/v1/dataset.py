import boto3
from flask import make_response, jsonify

from ....common.entities import Dataset
from ....common.utils.db_utils import db_session
from ....common.utils.exceptions import NotFoundHTTPException, CorporaException, ServerErrorHTTPException
from ....common.utils.s3_utils import generate_file_url, head_file


@db_session
def post_dataset_artifact(dataset_uuid: str, asset_uuid: str):
    s3 = boto3.client("s3")

    # retrieve the dataset
    dataset = Dataset.get(dataset_uuid)
    if not dataset:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}' not found.")

    # retrieve the artifact
    artifact = [artifact for artifact in dataset.artifacts if artifact.id == asset_uuid]
    if not artifact:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}/asset/{asset_uuid}' not found.")
    else:
        artifact = artifact[0]

    # Retrieve S3 metadata
    bucket_name, file_prefix = artifact.s3_uri[5:].split("/", 1)
    metadata = head_file(bucket_name, file_prefix, s3=s3)
    if not metadata:
        raise NotFoundHTTPException(f"'dataset/{dataset_uuid}/asset/{asset_uuid}' not found.")
    else:
        file_size = metadata["ContentLength"]

    # Generate presigned URL
    presigned_url = generate_file_url(bucket_name, file_prefix, s3=s3)
    if not presigned_url:
        raise ServerErrorHTTPException()

    return make_response(
        jsonify(
            dataset_id=dataset_uuid,
            file_name=artifact.filename,
            file_type=artifact.filetype,
            file_size=file_size,
            presigned_url=presigned_url,
        ),
        201,
    )
