import logging
import os
import re
import sys
from typing import Tuple, Optional
from urllib.parse import unquote_plus

from pythonjsonlogger import jsonlogger
from sqlalchemy.orm import Session

from backend.common.entities import Dataset, Collection
from backend.common.logging_config import DATETIME_FORMAT, LOG_FORMAT
from backend.common.upload import upload
from backend.common.utils.db_session import db_session_manager
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

        with db_session_manager() as session:
            collection_owner, dataset_id = get_dataset_info(session, parsed["collection_id"], parsed["dataset_id"])

            logger.info(dict(collection_owner=collection_owner, dataset_id=dataset_id))
            if not collection_owner:
                raise CorporaException(f"Collection {parsed['collection_id']} does not exist")
            elif parsed["username"] == "super":
                pass
            elif parsed["username"] != collection_owner:
                raise CorporaException(
                    f"user:{parsed['username']} does not have permission to modify datasets in collection "
                    f"{parsed['collection_id']}."
                )

            s3_uri = f"s3://{bucket}/{key}"
            upload(
                session,
                parsed["collection_id"],
                user=collection_owner,
                url=s3_uri,
                file_size=size,
                dataset_id=dataset_id,
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


def get_dataset_info(session: Session, collection_id: str, dataset_id: str) -> Tuple[Optional[str], Optional[str]]:
    if dataset := Dataset.get(session, dataset_id=dataset_id, collection_id=collection_id):
        return dataset.collection.owner, dataset.id
    else:
        if not Collection.get(session, collection_id):
            raise NonExistentCollectionException(f"No Collection with {collection_id=} found")
        raise NonExistentDatasetException(f"No Dataset with {dataset_id=} in Collection {collection_id}")
