import typing

from backend.corpora.common.utils.s3_utils import head_file, generate_file_url
from .entity import Entity
from ..corpora_orm import DbDatasetArtifact
import boto3


class Asset(Entity):
    table = DbDatasetArtifact
    s3 = boto3.client("s3")

    def __init__(self, db_object: DbDatasetArtifact):
        super().__init__(db_object)
        self.bucket_name, self.file_prefix = self.s3_uri[5:].split("/", 1)

    def generate_file_url(self, expiration: int = 3600) -> typing.Union[str, None]:
        return generate_file_url(self.bucket_name, self.file_prefix, expiration, self.s3)

    def get_s3_metadata(self) -> typing.Union[dict, None]:
        return head_file(self.bucket_name, self.file_prefix, s3=self.s3)

    def delete(self):
        """
        TODO: Must also delete s3 objects related to the asset.
        :return:
        """
        raise NotImplementedError()
