from unittest.mock import Mock

from backend.layers.thirdparty.s3_provider_mock import MockS3Provider
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProviderInterface
from backend.layers.thirdparty.uri_provider import FileInfo, UriProvider
from tests.unit.backend.layers.common.base_test import BaseTest


class BaseProcessingTest(BaseTest):
    def setUp(self):
        super().setUp()
        self.uri_provider = UriProvider()
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "local.h5ad"))
        self.s3_provider = MockS3Provider()
        self.schema_validator = Mock(spec=SchemaValidatorProviderInterface)
        self.schema_validator.validate_and_save_labels = Mock(return_value=(True, [], True))


mock_config_attr = {
    "collections_base_url": "http://collections",
    "dataset_assets_base_url": "http://domain",
    "upload_max_file_size_gb": 1,
    "schema_4_feature_flag": "True",
}


def mock_config_fn(name):
    return mock_config_attr[name]
