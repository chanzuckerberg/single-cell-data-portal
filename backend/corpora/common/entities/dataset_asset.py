import logging
import typing

import boto3
from botocore.exceptions import ClientError

from .entity import Entity
from ..corpora_orm import DbDatasetArtifact

logger = logging.getLogger(__name__)


class DatasetAsset(Entity):
    table = DbDatasetArtifact
    s3 = boto3.client("s3")

    def __init__(self, db_object: DbDatasetArtifact):
        super().__init__(db_object)

        self.bucket_name, self.file_prefix = self.s3_uri[5:].split("/", 1)

    def generate_file_url(self, expiration: int = 604800) -> typing.Union[str, None]:
        """
        Generate a presigned URL for a file for user download.
        :param expiration: Presigned URL expiration in seconds
        :return: Presigned URL to download the requested file
        """
        try:
            response = self.s3.generate_presigned_url(
                "get_object", Params={"Bucket": self.bucket_name, "Key": self.file_prefix}, ExpiresIn=expiration
            )
        except ClientError:
            logger.exception(f"Failed to generate presigned URL for '{self.file_prefix}'.")
            return None
        else:
            return response

    def get_file_size(self) -> typing.Union[int, None]:
        """
        Retrieves the asset content length from the S3 object.
        :return: The content length in bytes.
        """

        try:
            response = self.s3.head_object(Bucket=self.bucket_name, Key=self.file_prefix)
        except ClientError:
            logger.exception(f"Failed to retrieve meta data for '{self.file_prefix}'.")
            return None
        else:
            return response["ContentLength"]
