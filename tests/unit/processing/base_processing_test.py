from typing import List
from unittest.mock import Mock

from backend.layers.processing.downloader import Downloader
from backend.layers.thirdparty.s3_provider_interface import S3ProviderInterface
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProviderInterface
from backend.layers.thirdparty.uri_provider import FileInfo, UriProvider
from tests.unit.backend.layers.common.base_test import BaseTest


class MockS3Provider(S3ProviderInterface):
    """
    Simple S3 mock that mostly checks if paths are correct
    """

    def __init__(self) -> None:
        self.mock_s3_fs = []

    def upload_file(self, src_file: str, bucket_name: str, dst_file: str, extra_args: dict):
        url = f"s3://{bucket_name}/{dst_file}"
        self.mock_s3_fs.append(url)

    def upload_directory(self, src_dir: str, s3_uri: str):
        self.mock_s3_fs.append(s3_uri)

    def delete_files(self, bucket_name: str, object_keys: List[str]):
        pass

    def download_file(self, bucket_name: str, object_key: str, local_filename: str):
        pass

    def file_exists(self, bucket_name: str, object_key: str):
        url = f"s3://{bucket_name}/{object_key}"
        return self.uri_exists(url)

    def uri_exists(self, uri: str):
        return uri in self.mock_s3_fs

    def is_empty(self):
        return len(self.mock_s3_fs) == 0


class BaseProcessingTest(BaseTest):
    def setUp(self):
        super().setUp()
        self.uri_provider = UriProvider()
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "local.h5ad"))
        self.s3_provider = MockS3Provider()
        self.schema_validator = SchemaValidatorProviderInterface()
        self.schema_validator.validate_and_save_labels = Mock(return_value=(True, [], True))
        self.downloader = Downloader(self.business_logic)
        self.downloader.download_file = Mock()
