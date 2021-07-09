import os
from mock import patch
import json

from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestVersion(BaseAuthAPITest):
    @patch.dict(os.environ, {"COMMIT_SHA": "test"})
    def test_get(self):
        response = self.app.get("/dp/v1/deployed_version")
        self.assertEqual("test", json.loads(response.data)["Data Portal"])
