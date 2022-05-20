import logging
import os
import re
from typing import Tuple, Optional
from urllib.parse import unquote_plus

from sqlalchemy.orm import Session

from backend.corpora.common.entities import Dataset
from backend.corpora.common.upload import upload
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.exceptions import CorporaException

logger = logging.getLogger(__name__)
USERNAME_REGEX = r"(?P<username>[\w\-\|]+)"
UUID_REGEX = r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"
EXTENSION_REGEX = r"(?P<extension>h5ad)"
DATASET_ID_REGEX = f"(?P<dataset_uuid>{UUID_REGEX})"
COLLECTION_ID_REGEX = f"(?P<collection_uuid>{UUID_REGEX})"
REGEX = f"^{USERNAME_REGEX}/{COLLECTION_ID_REGEX}/((tag/(?P<tag>.*))|(id/{DATASET_ID_REGEX})).{EXTENSION_REGEX}$"


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
            raise CorporaException(f"Missing collection UUID, curator tag, and/or dataset UUID for {key=}")
        if parsed["tag"]:
            parsed["tag"] = f"{parsed['tag']}.{parsed['extension']}"
        logger.debug(parsed)

        with db_session_manager() as session:
            collection_owner, dataset_uuid = get_dataset_info(
                session, parsed["collection_uuid"], parsed["dataset_uuid"], parsed["tag"]
            )

            logger.debug(f"{collection_owner=}, {dataset_uuid=}")
            if not collection_owner:
                raise CorporaException(f"Collection {parsed['collection_uuid']} does not exist")
            elif parsed["username"] != collection_owner:
                raise CorporaException(
                    f"user:{parsed['username']} does not have permission to modify datasets in collection "
                    f"{parsed['collection_uuid']}."
                )

            s3_uri = f"s3://{bucket}/{key}"
            upload(
                session,
                parsed["collection_uuid"],
                user=collection_owner,
                url=s3_uri,
                file_size=size,
                file_extension=parsed["extension"],
                dataset_id=dataset_uuid,
                curator_tag=parsed["tag"],
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
    Parses the S3 object key to extract the collection UUID and curator tag, ignoring the REMOTE_DEV_PREFIX

    Example of key with only curator_tag:
    s3://<dataset submissions bucket>/<user_id>/<collection_id>/tag/<curator_tag>

    Example of key with dataset id:
    s3://<dataset submissions bucket>/<user_id>/<collection_id>/id/<dataset_id>

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
    session: Session, collection_uuid: str, dataset_uuid: str, incoming_curator_tag: str
) -> Tuple[Optional[str], Optional[str]]:
    if dataset_uuid:
        dataset = Dataset.get(session, dataset_uuid)
    else:
        dataset = Dataset.get_dataset_from_curator_tag(session, collection_uuid, incoming_curator_tag)
    if dataset:
        return dataset.collection.owner, dataset.id
    return None, None
