import typing

from backend.corpora.common.utils.s3_utils import head_file, generate_file_url
from .entity import Entity
from ..corpora_orm import DbDatasetArtifact


class DatasetAsset(Entity):
    table = DbDatasetArtifact

    def __init__(self, db_object: DbDatasetArtifact):
        super().__init__(db_object)
        self.bucket_name, self.file_prefix = self.s3_uri[5:].split("/", 1)

    def generate_file_url(self, expiration: int = 3600) -> typing.Union[str, None]:
        return generate_file_url(self.bucket_name, self.file_prefix, expiration)

    def get_file_size(self) -> typing.Union[dict, None]:
        metadata = head_file(self.bucket_name, self.file_prefix)
        return metadata["ContentLength"] if metadata else metadata
