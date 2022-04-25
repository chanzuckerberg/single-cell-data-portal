import boto3
import json
import time

from sqlalchemy.orm import Session

from .corpora_config import CorporaConfig
import os

from .corpora_orm import CollectionVisibility, ProcessingStatus
from .entities import Collection, Dataset
from .utils.exceptions import (
    MaxFileSizeExceededException,
    InvalidFileFormatException,
    NonExistentCollectionException,
    InvalidProcessingStateException,
    NonExistentDatasetException,
)
from .utils.math_utils import GB
from ..lambdas.api.v1.collection import _owner_or_allowed

_stepfunctions_client = None


def get_stepfunctions_client():
    global _stepfunctions_client
    if not _stepfunctions_client:
        _stepfunctions_client = boto3.client("stepfunctions", endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None)
    return _stepfunctions_client


def start_upload_sfn(collection_uuid, dataset_uuid, url):
    input_parameters = {"collection_uuid": collection_uuid, "url": url, "dataset_uuid": dataset_uuid}
    sfn_name = f"{dataset_uuid}_{int(time.time())}"
    response = get_stepfunctions_client().start_execution(
        stateMachineArn=CorporaConfig().upload_sfn_arn,
        name=sfn_name,
        input=json.dumps(input_parameters),
    )
    return response


def upload(
    db_session: Session,
    collection_uuid: str,
    user: str,
    url: str,
    file_size: int,
    file_extension: str,
    dataset_id: str = None,
) -> str:
    max_file_size_gb = CorporaConfig().upload_max_file_size_gb * GB
    if file_size is not None and file_size > max_file_size_gb:
        raise MaxFileSizeExceededException(f"{url} exceeds the maximum allowed file size of {max_file_size_gb} Gb")

    allowed_file_formats = CorporaConfig().upload_file_formats
    if file_extension not in allowed_file_formats:
        raise InvalidFileFormatException(f"{url} must be in the file format(s): {allowed_file_formats}")

    collection = Collection.get_collection(
        db_session,
        collection_uuid,
        visibility=CollectionVisibility.PRIVATE,  # Do not allow changes to public Collections
        owner=_owner_or_allowed(user),
    )
    if not collection:
        raise NonExistentCollectionException(f"Collection {collection_uuid} does not exist")

    if dataset_id:
        # Update dataset
        dataset = Dataset.get(db_session, dataset_id)
        if collection_uuid == dataset.collection_id:
            if dataset.processing_status.processing_status in [ProcessingStatus.SUCCESS, ProcessingStatus.FAILURE]:
                dataset.reprocess()
            else:
                raise InvalidProcessingStateException(
                    f"Unable to reprocess dataset {dataset_id}: {dataset.processing_status.processing_status=}"
                )
        else:
            raise NonExistentDatasetException(f"Dataset {dataset_id} does not exist")

    else:
        # Add new dataset
        dataset = Dataset.create(db_session, collection=collection)

    dataset.update(processing_status=dataset.new_processing_status())

    # Start processing link
    start_upload_sfn(collection_uuid, dataset.id, url)

    return dataset.id
