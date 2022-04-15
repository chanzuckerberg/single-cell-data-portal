import unittest
from unittest.mock import patch, Mock

from tests.unit.backend.corpora.api_server.base_api_test import BaseAPITest


class TestAuthToken(BaseAPITest):
    @patch("backend.corpora.lambdas.api.v1.curation.collection.dataset.sts_client")
    def test__generate_s3_credentials__OK(self, sts_client):
        pass


if __name__ == "__main__":
    unittest.main()
