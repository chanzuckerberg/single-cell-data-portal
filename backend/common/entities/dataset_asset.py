import logging
import typing
from os.path import basename, join
from urllib.parse import urlparse

from botocore.exceptions import ClientError

from backend.common.corpora_orm import DatasetArtifactFileType, DbDatasetArtifact
from backend.common.entities.entity import Entity
from backend.common.utils.s3_buckets import buckets


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
            response = buckets.portal_client.generate_presigned_url(
                "get_object", Params={"Bucket": self.bucket_name, "Key": self.key_name}, ExpiresIn=expiration
            )
        except ClientError:
            logging.exception(f"Failed to generate presigned URL for '{self.url}'.")
            return None
        else:
            return response

    def get_file_size(self) -> typing.Union[int, None]:
        """
        Retrieves the asset content length from the S3 object.
        :return: The content length in bytes.
        """

        try:
            response = buckets.portal_client.head_object(Bucket=self.bucket_name, Key=self.key_name)
        except ClientError:
            logging.exception(f"Failed to retrieve meta data for '{self.url}'.")
            return None
        else:
            return response["ContentLength"]

    def delete_from_s3(self):
        try:
            if self.key_name.endswith("/"):
                # This path should only be taken when deleting from the cellxgene bucket
                logging.info(f"Deleting all files in bucket {self.bucket_name} under {self.dataset_id}.")
                buckets.portal_resource.Bucket(self.bucket_name).objects.filter(Prefix=self.dataset_id).delete()
                # using dataset_id rather than the key_name because we also need to delete the genesets if they exist.
            else:
                logging.info(f"Deleting file {self.key_name} in bucket {self.bucket_name}.")
                buckets.portal_client.delete_object(Bucket=self.bucket_name, Key=self.key_name)
        except ClientError:
            logging.exception(f"Failed to delete artifact '{self.url}'.")

    def queue_s3_asset_for_deletion(self):
        """Deletes s3 assets after the database changes have been commited."""
        self.session.info["s3_deletion_list"].append(self.delete_from_s3)

    @classmethod
    def create(
        cls,
        session,
        dataset_id: str,
        filename: str,
        filetype: DatasetArtifactFileType,
        user_submitted: bool,
        s3_uri: str,
    ):

        db_object = cls.table(
            dataset_id=dataset_id,
            filename=filename,
            filetype=filetype,
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
        buckets.portal_client.upload_file(
            file_name,
            artifact_bucket,
            join(bucket_prefix, file_base),
            ExtraArgs={"ACL": "bucket-owner-full-control"},
        )
        return DatasetAsset.make_s3_uri(artifact_bucket, bucket_prefix, file_base)

    @classmethod
    def s3_uris_for_datasets(cls, session, dataset_ids, file_type=None, file_name=None) -> typing.Dict:
        """
        Returns a dictionary of dataset_id : s3_uri for a given list of dataset ids
        also filters on file_type if specified
        """
        s3_uris = session.query(cls.table.dataset_id, cls.table.s3_uri).filter(cls.table.dataset_id.in_(dataset_ids))
        if file_type:
            s3_uris = s3_uris.filter(cls.table.filetype == file_type)
        if file_name:
            s3_uris = s3_uris.filter(cls.table.filename == file_name)
        return {x[0]: x[1] for x in s3_uris.all()}
