import logging
import os
import re
from typing import Tuple, Optional
from urllib.parse import unquote_plus

from sqlalchemy.orm import Session

from backend.corpora.common.entities import Collection
from backend.corpora.common.upload import upload
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.exceptions import CorporaException

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
UUID_REGEX = "[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"


def dataset_submissions_handler(s3_event: dict, context) -> None:
    """
    Lambda function invoked when a dataset is uploaded to the dataset submissions S3 bucket
    :param s3_event: Lambda's event object
    :param context: Lambda's context object
    :return:
    """
    logger.debug(f"{s3_event=}")
    logger.debug(f"{os.environ.get('REMOTE_DEV_PREFIX', '')=}")

    # s3://<dataset submissions bucket>/<collection_id>/<curator_tag>

    for record in s3_event["Records"]:
        bucket, key, size = parse_s3_event_record(record)
        logger.debug(f"{bucket=}, {key=}, {size=}")

        collection_uuid, incoming_curator_tag = parse_key(key)

        if not collection_uuid or not incoming_curator_tag:
            raise CorporaException(f"Missing collection UUID and/or curator tag for {key=}")

        extension = get_extension(incoming_curator_tag)
        logger.debug(f"{collection_uuid=}, {incoming_curator_tag=}, {extension=}")

        with db_session_manager() as session:

            collection_owner, dataset_uuid = get_dataset_info(session, collection_uuid, incoming_curator_tag)

            logger.debug(f"{collection_owner=}, {dataset_uuid=}")
            if not collection_owner:
                raise CorporaException(f"Collection {collection_uuid} does not exist")

            s3_uri = f"s3://{bucket}/{key}"
            upload(
                session,
                collection_uuid,
                user=collection_owner,
                url=s3_uri,
                file_size=size,
                file_extension=extension,
                dataset_id=dataset_uuid,
                curator_tag=incoming_curator_tag,
                is_api_call=False,
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


def parse_key(key: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parses the S3 object key to extract the collection UUID and curator tag, ignoring the REMOTE_DEV_PREFIX
    :param key:
    :return:
    """
    rdev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "").strip("/")
    if rdev_prefix:
        key = key.replace(f"{rdev_prefix}/", "")

    matched = re.match(f"^({UUID_REGEX})/(.*)$", key)
    if matched:
        collection_uuid, curator_tag = matched.groups()
        return collection_uuid.lower(), curator_tag
    else:
        return None, None


def get_extension(path: str) -> str:
    extension = os.path.splitext(path)[1]
    return extension.replace(".", "")


def get_dataset_info(
    session: Session, collection_uuid: str, incoming_curator_tag: str
) -> Tuple[Optional[str], Optional[str]]:
    collection = Collection.get_collection(session=session, collection_uuid=collection_uuid)

    dataset_uuid = None

    if not collection:
        return None, dataset_uuid

    for dataset in collection.datasets:
        if dataset.curator_tag == incoming_curator_tag:
            dataset_uuid = dataset.id
            break

    return collection.owner, dataset_uuid
