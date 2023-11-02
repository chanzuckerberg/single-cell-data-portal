import logging
from datetime import datetime
from os.path import basename, join
from typing import Callable, List, Optional

from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    DatasetConversionStatus,
    DatasetStatusGeneric,
    DatasetStatusKey,
    DatasetVersion,
    DatasetVersionId,
)
from backend.layers.processing.exceptions import ConversionFailed
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
        dataset_id: DatasetVersionId,
        status_key: DatasetStatusKey,
        status_value: DatasetStatusGeneric,
        validation_errors: Optional[List[str]] = None,
    ):
        validation_message = "\n".join(validation_errors) if validation_errors is not None else None
        self.business_logic.update_dataset_version_status(dataset_id, status_key, status_value, validation_message)
        self.logger.info(
            "Updating processing status",
            extra=dict(validation_message=validation_message, status_key=status_key, status_value=status_value),
        )

    def download_from_s3(self, bucket_name: str, object_key: str, local_filename: str):
        self.s3_provider.download_file(bucket_name, object_key, local_filename)

    @staticmethod
    def make_s3_uri(artifact_bucket, key_prefix, file_name):
        return join("s3://", artifact_bucket, key_prefix, file_name)

    def upload_artifact(
        self,
        file_name: str,
        key_prefix: str,
        artifact_bucket: str,
    ) -> str:
        file_base = basename(file_name)
        self.s3_provider.upload_file(
            file_name,
            artifact_bucket,
            join(key_prefix, file_base),
            extra_args={"ACL": "bucket-owner-full-control"},
        )
        return self.make_s3_uri(artifact_bucket, key_prefix, file_base)

    @logit
    def create_artifact(
        self,
        file_name: str,
        artifact_type: str,
        key_prefix: str,
        dataset_id: DatasetVersionId,
        artifact_bucket: str,
        processing_status_key: DatasetStatusKey,
        datasets_bucket: Optional[str] = None,  # If provided, dataset will be uploaded to this bucket for public access
    ):
        self.update_processing_status(dataset_id, processing_status_key, DatasetConversionStatus.UPLOADING)
        try:
            s3_uri = self.upload_artifact(file_name, key_prefix, artifact_bucket)
            self.logger.info(f"Uploaded [{dataset_id}/{file_name}] to {s3_uri}")
            self.business_logic.add_dataset_artifact(dataset_id, artifact_type, s3_uri)
            self.logger.info(f"Updated database with {artifact_type}.")
            if datasets_bucket:
                key = ".".join((key_prefix, artifact_type))
                self.s3_provider.upload_file(
                    file_name, datasets_bucket, key, extra_args={"ACL": "bucket-owner-full-control"}
                )
                datasets_s3_uri = self.make_s3_uri(datasets_bucket, key_prefix, key)
                self.logger.info(f"Uploaded {dataset_id}.{artifact_type} to {datasets_s3_uri}")
            self.update_processing_status(dataset_id, processing_status_key, DatasetConversionStatus.UPLOADED)
        except Exception:
            raise ConversionFailed(processing_status_key) from None

    def convert_file(
        self,
        converter: Callable,
        local_filename: str,
        error_message: str,
        dataset_id: DatasetVersionId,
        processing_status_key: DatasetStatusKey,
    ) -> str:
        self.logger.info(f"Converting {local_filename}")
        start = datetime.now()
        try:
            self.update_processing_status(dataset_id, processing_status_key, DatasetConversionStatus.CONVERTING)
            file_dir = converter(local_filename)
            self.update_processing_status(dataset_id, processing_status_key, DatasetConversionStatus.CONVERTED)
            self.logger.info(f"Finished converting {converter} in {datetime.now() - start}")
        except Exception:
            raise ConversionFailed(processing_status_key) from None
        return file_dir

    def get_key_prefix(self, identifier: str) -> str:
        import os

        remote_dev_prefix = os.environ.get("REMOTE_DEV_PREFIX", "")
        if remote_dev_prefix:
            return join(remote_dev_prefix, identifier).strip("/")
        else:
            return identifier

    def check_dataset_is_latest_schema_version(self, dataset: DatasetVersion) -> bool:
        return (
            hasattr(dataset, "metadata")
            and hasattr(dataset.metadata, "schema_version")
            and dataset.metadata.schema_version == self.schema_version
        )
