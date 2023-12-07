import logging
import os
import re
import sys
from typing import Optional, Tuple
from urllib.parse import unquote_plus

from pythonjsonlogger import jsonlogger

from backend.common.logging_config import DATETIME_FORMAT, LOG_FORMAT
from backend.common.utils.exceptions import (
    CorporaException,
    NonExistentCollectionException,
    NonExistentDatasetException,
)
from backend.common.utils.regex import COLLECTION_ID_REGEX, DATASET_ID_REGEX, USERNAME_REGEX
from backend.layers.business.business import BusinessLogic
from backend.layers.business.exceptions import CollectionNotFoundException, DatasetNotFoundException
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.batch_job_provider import BatchJobProvider
from backend.layers.thirdparty.crossref_provider import CrossrefProvider
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.step_function_provider import StepFunctionProvider
from backend.layers.thirdparty.uri_provider import UriProvider

log_handler = logging.StreamHandler(stream=sys.stdout)
formatter = jsonlogger.JsonFormatter(LOG_FORMAT, DATETIME_FORMAT)
log_handler.setFormatter(formatter)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.handlers = [log_handler]

REGEX = f"^{USERNAME_REGEX}/{COLLECTION_ID_REGEX}/{DATASET_ID_REGEX}$"

_business_logic = None


def get_business_logic():
    global _business_logic
    if not _business_logic:
        database_provider = DatabaseProvider()
        uri_provider = UriProvider()
        step_function_provider = StepFunctionProvider()
        s3_provider = S3Provider()
        crossref_provider = CrossrefProvider()
        batch_job_provider = BatchJobProvider()
        _business_logic = BusinessLogic(
            database_provider, batch_job_provider, crossref_provider, step_function_provider, s3_provider, uri_provider
        )
    return _business_logic


def dataset_submissions_handler(s3_event: dict, unused_context) -> None:
    """
    Lambda function invoked when a dataset is uploaded to the dataset submissions S3 bucket
    :param s3_event: Lambda's event object
    :param unused_context: Lambda's context object
    :return:
    """
    logger.info(dict(message="s3_event", **s3_event))
    logger.debug(dict(REMOTE_DEV_PREFIX=os.environ.get("REMOTE_DEV_PREFIX", "")))

    for record in s3_event["Records"]:
        bucket, key, size = parse_s3_event_record(record)
        logger.debug(f"{bucket=}, {key=}, {size=}")

        parsed = parse_key(key)
        if not parsed:
            raise CorporaException(f"Missing Collection ID and/or Dataset ID for {key=}")
        logger.debug(parsed)

        collection_id = parsed["collection_id"]
        dataset_id = parsed["dataset_id"]

        business_logic = get_business_logic()
        try:
            collection_version, dataset_version = business_logic._get_collection_and_dataset(collection_id, dataset_id)
        except CollectionNotFoundException:
            raise NonExistentCollectionException(f"Collection {parsed['collection_id']} does not exist") from None
        except DatasetNotFoundException:
            raise NonExistentDatasetException(f"No Dataset with {dataset_id=} in Collection {collection_id}") from None

        collection_owner = collection_version.owner

        logger.info(dict(collection_owner=collection_owner, dataset_id=dataset_id))
        if parsed["username"] == "super":
            pass
        elif parsed["username"] != collection_owner:
            raise CorporaException(
                f"user:{parsed['username']} does not have permission to modify datasets in collection "
                f"{parsed['collection_id']}."
            )

        s3_uri = f"s3://{bucket}/{key}"

        get_business_logic().ingest_dataset(
            collection_version_id=collection_version.version_id,
            url=s3_uri,
            file_size=size,
            current_dataset_version_id=dataset_version.version_id,
        )


def parse_s3_event_record(s3_event_record: dict) -> Tuple[str, str, int]:
    """
    Parses the S3 event record and returns the bucket name, object key and object size
    :param s3_event_record:
    :return:
    """
    bucket = s3_event_record["s3"]["bucket"]["name"]
    key = unquote_plus(s3_event_record["s3"]["object"]["key"], encoding="utf-8")
    size = s3_event_record["s3"]["object"]["size"]
    return bucket, key, size


def parse_key(key: str) -> Optional[dict]:
    """
    Parses the S3 object key to extract the Collection ID and Dataset ID, ignoring the REMOTE_DEV_PREFIX

    Example of key with dataset id:
    s3://<dataset submissions bucket>/<user_id>/<collection_id>/<dataset_id>

    :param key:
    :return:
    """
    rdev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "").strip("/")
    if rdev_prefix:
        key = key.replace(f"{rdev_prefix}/", "")

    matched = re.match(REGEX, key)
    if matched:
        return matched.groupdict()
