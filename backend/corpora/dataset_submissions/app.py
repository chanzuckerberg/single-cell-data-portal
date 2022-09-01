import logging
import os
import re
from typing import Tuple, Optional
from urllib.parse import unquote_plus

from sqlalchemy.orm import Session

from backend.corpora.common.entities import Dataset, Collection
from backend.corpora.common.upload import upload
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.exceptions import CorporaException
from backend.corpora.common.utils.regex import (
    USERNAME_REGEX,
    COLLECTION_ID_REGEX,
    DATASET_ID_REGEX,
    CURATOR_TAG_REGEX,
)

logger = logging.getLogger(__name__)
REGEX = f"^{USERNAME_REGEX}/{COLLECTION_ID_REGEX}/({DATASET_ID_REGEX}|{CURATOR_TAG_REGEX})$"


def dataset_submissions_handler(s3_event: dict, unused_context) -> None:
    """
    Lambda function invoked when a dataset is uploaded to the dataset submissions S3 bucket
    :param s3_event: Lambda's event object
    :param unused_context: Lambda's context object
    :return:
    """
    logger.debug(f"{s3_event=}")
    logger.debug(f"{os.environ.get('REMOTE_DEV_PREFIX', '')=}")

    for record in s3_event["Records"]:
        bucket, key, size = parse_s3_event_record(record)
        logger.debug(f"{bucket=}, {key=}, {size=}")

        parsed = parse_key(key)
        if not parsed:
            raise CorporaException(f"Missing collection ID, curator tag, and/or dataset ID for {key=}")
        logger.debug(parsed)

        with db_session_manager() as session:
            collection_owner, dataset_id = get_dataset_info(
                session, parsed["collection_id"], parsed["dataset_id"], parsed.get("tag")
            )

            logger.info(f"{collection_owner=}, {dataset_id=}")
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
                curator_tag=parsed.get("tag"),
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
    Parses the S3 object key to extract the collection ID and curator tag, ignoring the REMOTE_DEV_PREFIX

    Example of key with only curator_tag:
    s3://<dataset submissions bucket>/<user_id>/<collection_id>/<curator_tag>

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


def get_dataset_info(
    session: Session, collection_id: str, dataset_id: str, incoming_curator_tag: str
) -> Tuple[Optional[str], Optional[str]]:
    if dataset_id:  # If a dataset uuid was provided
        dataset = Dataset.get(session, dataset_id)
    elif incoming_curator_tag:  # if incoming_curator_tag
        dataset = Dataset.get_dataset_from_curator_tag(session, collection_id, incoming_curator_tag)
        if not dataset:  # New dataset
            collection = Collection.get_collection(session, collection_id)
            if collection:
                return collection.owner, None
    else:
        raise CorporaException("No dataset identifier provided")
    if dataset:
        return dataset.collection.owner, dataset.id
    return None, None
