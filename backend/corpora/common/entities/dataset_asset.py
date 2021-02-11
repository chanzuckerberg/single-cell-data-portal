import logging
import os
import typing
from urllib.parse import urlparse

import boto3
from botocore.exceptions import ClientError

from .entity import Entity
from ..corpora_orm import DbDatasetArtifact, DatasetArtifactType, DatasetArtifactFileType

logger = logging.getLogger(__name__)


class DatasetAsset(Entity):
    table = DbDatasetArtifact

    _s3 = None

    @classmethod
    def s3_client(cls):
        if not cls._s3:
            cls._s3 = boto3.client(
                "s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"), config=boto3.session.Config(signature_version="s3v4")
            )
        return cls._s3

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
            response = self.s3_client().generate_presigned_url(
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
            response = self.s3_client().head_object(Bucket=self.bucket_name, Key=self.key_name)
        except ClientError:
            logger.exception(f"Failed to retrieve meta data for '{self.url}'.")
            return None
        else:
            return response["ContentLength"]

    def delete_from_s3(self):
        try:
            self.s3_client().delete_object(Bucket=self.bucket_name, Key=self.key_name)
        except ClientError:
            logger.exception(f"Failed to delete artifact '{self.url}'.")
            return None

    @classmethod
    def create(
        cls,
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
        cls.db.session.add(db_object)
        cls.db.commit()
        return cls(db_object)
