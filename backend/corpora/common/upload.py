import boto3
import json
import time

from sqlalchemy.orm import Session

from backend.corpora.common.corpora_config import CorporaConfig
import os

from backend.corpora.common.corpora_orm import CollectionVisibility, ProcessingStatus, ValidationStatus
from backend.corpora.common.entities import Collection, Dataset
from backend.corpora.common.utils.authorization_checks import owner_or_allowed
from backend.corpora.common.utils.exceptions import (
    MaxFileSizeExceededException,
    NonExistentCollectionException,
    InvalidProcessingStateException,
    NonExistentDatasetException,
)
from backend.corpora.common.utils.math_utils import GB

_stepfunctions_client = None


def get_stepfunctions_client():
    global _stepfunctions_client
    if not _stepfunctions_client:
        _stepfunctions_client = boto3.client("stepfunctions", endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None)
    return _stepfunctions_client


def start_upload_sfn(collection_id, dataset_id, url):
    input_parameters = {
        "collection_id": collection_id,
        "url": url,
        "dataset_id": dataset_id,
    }
    sfn_name = f"{dataset_id}_{int(time.time())}"
    response = get_stepfunctions_client().start_execution(
        stateMachineArn=CorporaConfig().upload_sfn_arn,
        name=sfn_name,
        input=json.dumps(input_parameters),
    )
    return response


def upload(
    db_session: Session,
    collection_id: str,
    url: str,
    file_size: int,
    user: str,
    scope: str = None,
    dataset_id: str = None,
) -> str:
    # Check if a dataset already exists
    if dataset_id:
        dataset = Dataset.get(db_session, dataset_id, collection_id=collection_id)
        if not dataset:
            raise NonExistentDatasetException(f"Dataset {dataset_id} does not exist")
    else:
        dataset = None

    # Check if file size is within max size limit
    max_file_size_gb = CorporaConfig().upload_max_file_size_gb
    if file_size is not None and file_size > max_file_size_gb * GB:
        error_message = f"The file exceeds the maximum allowed size of {max_file_size_gb} GB"
        if dataset:
            dataset.update(
                processing_status={
                    "validation_status": ValidationStatus.INVALID,
                    "validation_message": error_message,
                }
            )
        raise MaxFileSizeExceededException(error_message)

    # Check if datasets can be added to the collection
    collection = Collection.get_collection(
        db_session,
        collection_id,
        visibility=CollectionVisibility.PRIVATE,  # Do not allow changes to public Collections
        owner=owner_or_allowed(user, scope) if scope else user,
    )
    if not collection:
        raise NonExistentCollectionException(f"Collection {collection_id} does not exist")

    if dataset:
        # Update dataset
        if dataset.processing_status.processing_status not in [
            ProcessingStatus.SUCCESS,
            ProcessingStatus.FAILURE,
            ProcessingStatus.INITIALIZED,
        ]:
            raise InvalidProcessingStateException(
                f"Unable to reprocess dataset {dataset_id}: {dataset.processing_status.processing_status=}"
            )
        else:
            dataset.reprocess()

    else:
        # Add new dataset
        dataset = Dataset.create(db_session, collection=collection)

    dataset.update(processing_status=dataset.new_processing_status())

    # Start processing link
    start_upload_sfn(collection_id, dataset.id, url)

    return dataset.id
