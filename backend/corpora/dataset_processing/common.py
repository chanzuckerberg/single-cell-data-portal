import os
import typing
from datetime import datetime
from os.path import join

from backend.corpora.common.corpora_orm import DatasetArtifactFileType, ConversionStatus
from backend.corpora.common.entities import DatasetAsset, Dataset
from backend.corpora.common.utils.db_helpers import processing_status_updater
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.s3_buckets import buckets
from backend.corpora.dataset_processing.exceptions import ProcessingCancelled, ConversionFailed
from backend.corpora.dataset_processing.logger import logger


# This is unfortunate, but this information doesn't appear to live anywhere
# accessible to the uploader
DEPLOYMENT_STAGE_TO_URL = {
    "staging": "https://cellxgene.staging.single-cell.czi.technology/e",
    "prod": "https://cellxgene.cziscience.com/e",
    "rdev": f"https:/{os.environ.get('REMOTE_DEV_PREFIX')}-explorer.rdev.single-cell.czi.technology/e",
    "dev": "https://cellxgene.dev.single-cell.czi.technology/e",
    "test": "http://frontend.corporanet.local:3000",
}


def create_artifact(
    file_name: str,
    artifact_type: DatasetArtifactFileType,
    bucket_prefix: str,
    dataset_id: str,
    artifact_bucket: str,
    processing_status_type: str = None,
):
    if processing_status_type:
        update_db(
            dataset_id,
            processing_status={processing_status_type: ConversionStatus.UPLOADING},
        )
    logger.info(f"Uploading [{dataset_id}/{file_name}] to S3 bucket: [{artifact_bucket}].")
    try:
        s3_uri = DatasetAsset.upload(file_name, bucket_prefix, artifact_bucket)
        with db_session_manager() as session:
            logger.info(f"Updating database with  {artifact_type}.")
            DatasetAsset.create(
                session,
                dataset_id=dataset_id,
                filename=file_name,
                filetype=artifact_type,
                user_submitted=True,
                s3_uri=s3_uri,
            )
        if processing_status_type:
            update_db(
                dataset_id,
                processing_status={processing_status_type: ConversionStatus.UPLOADED},
            )

    except Exception as e:
        logger.error(e)
        e.args = {processing_status_type: ConversionStatus.FAILED}
        raise e


def update_db(dataset_id, metadata: dict = None, processing_status: dict = None):
    with db_session_manager() as session:
        dataset = Dataset.get(session, dataset_id, include_tombstones=True)
        if dataset.tombstone:
            raise ProcessingCancelled

        if metadata:
            logger.debug("Updating metadata.")
            dataset.update(**metadata)

        if processing_status:
            logger.debug(f"updating processing_status.{processing_status}")
            processing_status_updater(session, dataset.processing_status.id, processing_status)


def download_from_s3(
    bucket_name: str, object_key: str, local_filename: str, callback: typing.Optional[typing.Callable[int]] = None
):
    logger.info(f"Downloading file {local_filename} from bucket {bucket_name} with object key {object_key}")
    buckets.portal_client.download_file(bucket_name, object_key, local_filename, callback=None)


def convert_file(
    converter: typing.Callable,
    local_filename: str,
    error_message: str,
    dataset_id: str,
    processing_status_type: str,
) -> str:
    logger.info(f"Converting {local_filename}")
    start = datetime.now()
    try:
        update_db(
            dataset_id,
            processing_status={processing_status_type: ConversionStatus.CONVERTING},
        )
        file_dir = converter(local_filename)
        update_db(
            dataset_id,
            processing_status={processing_status_type: ConversionStatus.CONVERTED},
        )
        logger.info(f"Finished converting {converter} in {datetime.now()- start}")
    except Exception as e:
        logger.exception(f"{error_message}: {e}")
        status = {processing_status_type: ConversionStatus.FAILED}
        raise ConversionFailed(status)
    return file_dir


def get_bucket_prefix(identifier):
    remote_dev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "")
    if remote_dev_prefix:
        return join(remote_dev_prefix, identifier).strip("/")
    else:
        return identifier
