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
        self.schema_validator.validate_anndata = Mock(return_value=(True, [], True))
        self.schema_validator.validate_atac = Mock(return_value=([], "fragment.tsv.bgz.tbi", "fragment.tsv.bgz"))
