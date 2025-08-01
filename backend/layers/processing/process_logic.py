import logging
import os
from datetime import datetime
from os.path import basename
from typing import Callable, List, Optional

from backend.common.utils.dl_sources.uri import DownloadFailed
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    ARTIFACT_TO_EXTENSION,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetStatusGeneric,
    DatasetStatusKey,
    DatasetVersion,
    DatasetVersionId,
)
from backend.layers.processing.exceptions import ConversionFailed, UploadFailed
from backend.layers.processing.logger import logit
from backend.layers.thirdparty.s3_provider import S3ProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface


class ProcessingLogic:  # TODO: ProcessingLogicBase
    """
    Base class that contains all the processing logic methods
    """

    business_logic: BusinessLogicInterface
    uri_provider: UriProviderInterface
    s3_provider: S3ProviderInterface
    logger: logging.Logger
    schema_version: str

    def __init__(self) -> None:
        self.logger = logging.getLogger("processing")

    def update_processing_status(
        self,
        dataset_version_id: DatasetVersionId,
        status_key: DatasetStatusKey,
        status_value: DatasetStatusGeneric,
        validation_errors: Optional[List[str]] = None,
    ):
        validation_message = "\n".join(validation_errors) if validation_errors is not None else None
        self.business_logic.update_dataset_version_status(
            dataset_version_id, status_key, status_value, validation_message
        )
        self.logger.info(
            "Updating processing status",
            extra=dict(
                validation_message=validation_message,
                status_key=status_key,
                status_value=status_value,
                dataset_version_id=dataset_version_id.id,
            ),
        )

    @logit
    def download_from_source_uri(self, source_uri: str, local_path: str) -> str:
        """Given a source URI, download it to local_path.
        Handles fixing the url so it downloads directly.
        """
        file_url = self.uri_provider.parse(source_uri)
        if not file_url:
            raise ValueError(f"Malformed source URI: {source_uri}")
        try:
            file_url.download(local_path)
        except DownloadFailed as e:
            raise UploadFailed(f"Failed to download file from source URI: {source_uri}") from e
        return local_path

    def download_from_s3(self, bucket_name: str, object_key: str, local_filename: str):
        self.s3_provider.download_file(bucket_name, object_key, local_filename)

    def upload_artifact(self, file_name: str, key: str, bucket_name: str) -> str:
        self.s3_provider.upload_file(
            file_name,
            bucket_name,
            key,
            extra_args={"ACL": "bucket-owner-full-control"},
        )
        return "/".join(["s3:/", bucket_name, key])

    @logit
    def create_artifact(
        self,
        file_name: str,
        artifact_type: DatasetArtifactType,
        key_prefix: str,
        dataset_version_id: DatasetVersionId,
        processing_status_key: DatasetStatusKey,
        artifact_bucket: str,  # If provided, dataset will be uploaded to this bucket for future migrations
        datasets_bucket: Optional[str] = None,  # If provided, dataset will be uploaded to this bucket for public access
    ):
        self.update_processing_status(dataset_version_id, processing_status_key, DatasetConversionStatus.UPLOADING)
        try:
            key = "/".join([key_prefix, basename(file_name)])
            s3_uri = self.upload_artifact(file_name, key, artifact_bucket)
            self.logger.info(f"Uploaded [{dataset_version_id}/{file_name}] to {s3_uri}")
            self.business_logic.add_dataset_artifact(dataset_version_id, artifact_type, s3_uri)
            self.logger.info(f"Updated database with {artifact_type}.")
            if datasets_bucket:
                key = ".".join([key_prefix, ARTIFACT_TO_EXTENSION[artifact_type]])
                s3_uri = self.upload_artifact(file_name, key, datasets_bucket)
                self.logger.info(f"Uploaded [{dataset_version_id}/{file_name}] to {s3_uri}")
            self.update_processing_status(dataset_version_id, processing_status_key, DatasetConversionStatus.UPLOADED)
        except Exception as e:
            self.logger.error(e)
            raise ConversionFailed(processing_status_key) from None

    def convert_file(
        self,
        converter: Callable,
        local_filename: str,
        dataset_version_id: DatasetVersionId,
        fragment_artifact_id: Optional[str],
        processing_status_key: DatasetStatusKey,
    ) -> str:
        self.logger.info(f"Converting {local_filename}")
        start = datetime.now()
        try:
            self.update_processing_status(dataset_version_id, processing_status_key, DatasetConversionStatus.CONVERTING)
            file_dir = converter(local_filename, dataset_version_id, fragment_artifact_id)
            self.update_processing_status(dataset_version_id, processing_status_key, DatasetConversionStatus.CONVERTED)
            self.logger.info(f"Finished converting {converter} in {datetime.now() - start}")
        except Exception:
            raise ConversionFailed(processing_status_key) from None
        return file_dir

    def get_key_prefix(self, identifier: str) -> str:

        remote_dev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "")
        if remote_dev_prefix:
            return "/".join([remote_dev_prefix, identifier]).strip("/")
        else:
            return identifier

    def check_dataset_is_latest_schema_version(self, dataset: DatasetVersion) -> bool:
        return (
            hasattr(dataset, "metadata")
            and hasattr(dataset.metadata, "schema_version")
            and dataset.metadata.schema_version == self.schema_version
        )
