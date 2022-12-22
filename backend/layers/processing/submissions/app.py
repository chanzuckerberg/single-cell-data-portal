import logging
import os
import re
import sys
from typing import Tuple, Optional
from urllib.parse import unquote_plus
from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.common.entities import CollectionVersionId, DatasetVersionId

from pythonjsonlogger import jsonlogger

from backend.common.logging_config import DATETIME_FORMAT, LOG_FORMAT
from backend.common.utils.exceptions import (
    CorporaException,
    NonExistentCollectionException,
    NonExistentDatasetException,
)
from backend.common.utils.regex import USERNAME_REGEX, COLLECTION_ID_REGEX, DATASET_ID_REGEX

log_handler = logging.StreamHandler(stream=sys.stdout)
formatter = jsonlogger.JsonFormatter(LOG_FORMAT, DATETIME_FORMAT)
log_handler.setFormatter(formatter)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.handlers = [log_handler]

REGEX = f"^{USERNAME_REGEX}/{COLLECTION_ID_REGEX}/{DATASET_ID_REGEX}$"

_business_logic: BusinessLogic


def get_business_logic():
    global _business_logic
    if _business_logic is None:
        database_provider = DatabaseProvider()
        _business_logic = BusinessLogic(database_provider, None, None, None, None)
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

        collection_version_id = CollectionVersionId(parsed["collection_id"])
        dataset_version_id = DatasetVersionId(parsed["dataset_id"])

        version = get_business_logic().get_collection_version(collection_version_id)
        if version is None:
            raise NonExistentCollectionException(f"Collection {parsed['collection_id']} does not exist")

        if dataset_version_id not in [d.version_id for d in version.datasets]:
            raise NonExistentDatasetException(
                f"No Dataset with {dataset_version_id=} in Collection {collection_version_id}"
            )

        collection_owner = version.owner

        logger.info(dict(collection_owner=collection_owner, dataset_id=dataset_version_id))
        if parsed["username"] == "super":
            pass
        elif parsed["username"] != collection_owner:
            raise CorporaException(
                f"user:{parsed['username']} does not have permission to modify datasets in collection "
                f"{parsed['collection_id']}."
            )

        s3_uri = f"s3://{bucket}/{key}"

        get_business_logic().ingest_dataset(
            collection_version_id=collection_version_id,
            url=s3_uri,
            file_size=size,
            existing_dataset_version_id=dataset_version_id,
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
