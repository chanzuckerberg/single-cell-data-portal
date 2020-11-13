import logging
import typing
import os

import boto3
from botocore.exceptions import ClientError

from .entity import Entity
from ..corpora_orm import DbDatasetArtifact
from urllib.parse import urlparse

logger = logging.getLogger(__name__)


class DatasetAsset(Entity):
    table = DbDatasetArtifact
    s3 = boto3.client(
        "s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"), config=boto3.session.Config(signature_version="s3v4")
    )

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
            response = self.s3.generate_presigned_url(
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
            response = self.s3.head_object(Bucket=self.bucket_name, Key=self.key_name)
        except ClientError:
            logger.exception(f"Failed to retrieve meta data for '{self.url}'.")
            return None
        else:
            return response["ContentLength"]
