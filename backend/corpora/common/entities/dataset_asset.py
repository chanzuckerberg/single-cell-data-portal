import logging
from os.path import basename, join

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

    _s3_client = None
    _s3_resource = None

    @classmethod
    def s3_client(cls):
        if not cls._s3_client:
            cls._s3_client = boto3.client(
                "s3",
                endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None,
                config=boto3.session.Config(signature_version="s3v4"),
            )
        return cls._s3_client

    @classmethod
    def s3_resource(cls):
        if not cls._s3_resource:
            cls._s3_resource = boto3.resource(
                "s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"), config=boto3.session.Config(signature_version="s3v4")
            )
        return cls._s3_resource

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
            if self.key_name.endswith("/"):
                # This path should only be taken when deleting from the cellxgene bucket
                logger.info(f"Deleting all files in bucket {self.bucket_name} under {self.dataset_id}.")
                self.s3_resource().Bucket(self.bucket_name).objects.filter(Prefix=self.dataset_id).delete()
                # using dataset_id rather than the key_name because we also need to delete the genesets if they exist.
            else:
                logger.info(f"Deleting file {self.key_name} in bucket {self.bucket_name}.")
                self.s3_client().delete_object(Bucket=self.bucket_name, Key=self.key_name)
        except ClientError:
            logger.exception(f"Failed to delete artifact '{self.url}'.")

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

    def get_bucket_path(self):
        return "/".join(self.s3_uri.split("/")[-2:])

    @staticmethod
    def make_s3_uri(artifact_bucket, bucket_prefix, file_name):
        return join("s3://", artifact_bucket, bucket_prefix, file_name)

    @classmethod
    def upload(
        cls,
        file_name: str,
        bucket_prefix: str,
        artifact_bucket: str,
    ) -> str:
        file_base = basename(file_name)
        cls.s3_client().upload_file(
            file_name,
            artifact_bucket,
            join(bucket_prefix, file_base),
            ExtraArgs={"ACL": "bucket-owner-full-control"},
        )
        return DatasetAsset.make_s3_uri(artifact_bucket, bucket_prefix, file_base)
