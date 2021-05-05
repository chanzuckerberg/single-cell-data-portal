from os.path import basename, join

import logging
import typing
from urllib.parse import urlparse

from botocore.exceptions import ClientError

from .entity import Entity
from ..corpora_orm import DbDatasetArtifact, DatasetArtifactType, DatasetArtifactFileType
from ..utils.s3_buckets import s3_client

logger = logging.getLogger(__name__)


class DatasetAsset(Entity):
    table = DbDatasetArtifact

    def __init__(self, db_object: DbDatasetArtifact):
        super().__init__(db_object)

        self.url = urlparse(self.s3_uri)
        self.bucket_name = self.url.netloc
        self.key_name = self.url.path[1:]

    def generate_file_url(self, expiration: int = 604800) -> typing.Union[str, None]:
        """
        Generate a presigned URL for a file for user download.
        :param expiration: Presigned URL expiration in seconds. The default is 1 week.
        :return: Presigned URL to download the requested file
        """
        try:
            response = s3_client.generate_presigned_url(
                "get_object", Params={"Bucket": self.bucket_name, "Key": self.key_name}, ExpiresIn=expiration
            )
        except ClientError:
            logger.exception(f"Failed to generate presigned URL for '{self.url}'.")
            return None
        else:
            return response

    def get_file_size(self) -> typing.Union[int, None]:
        """
        Retrieves the asset content length from the S3 object.
        :return: The content length in bytes.
        """

        try:
            response = s3_client.head_object(Bucket=self.bucket_name, Key=self.key_name)
        except ClientError:
            logger.exception(f"Failed to retrieve meta data for '{self.url}'.")
            return None
        else:
            return response["ContentLength"]

    def delete_from_s3(self):
        try:
            s3_client.delete_object(Bucket=self.bucket_name, Key=self.key_name)
        except ClientError:
            logger.exception(f"Failed to delete artifact '{self.url}'.")
            return None

    @classmethod
    def create(
        cls,
        session,
        dataset_id: str,
        filename: str,
        filetype: DatasetArtifactFileType,
        type_enum: DatasetArtifactType,
        user_submitted: bool,
        s3_uri: str,
    ):

        db_object = cls.table(
            dataset_id=dataset_id,
            filename=filename,
            filetype=filetype,
            type=type_enum,
            user_submitted=user_submitted,
            s3_uri=s3_uri,
        )
        session.add(db_object)
        session.commit()
        return cls(db_object)

    @staticmethod
    def make_s3_uri(artifact_bucket, bucket_prefix, file_name):
        return join("s3://", artifact_bucket, bucket_prefix, file_name)

    @staticmethod
    def upload(
        file_name: str,
        bucket_prefix: str,
        artifact_bucket: str,
    ) -> str:
        file_base = basename(file_name)
        logger.info(f"Uploading to [{file_base}] to S3 bucket: [{artifact_bucket}].")
        s3_client.upload_file(
            file_name,
            artifact_bucket,
            join(bucket_prefix, file_base),
            ExtraArgs={"ACL": "bucket-owner-full-control"},
        )
        return DatasetAsset.make_s3_uri(artifact_bucket, bucket_prefix, file_base)
