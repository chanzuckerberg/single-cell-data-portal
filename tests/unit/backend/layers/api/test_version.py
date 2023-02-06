import json
import os
from unittest.mock import patch

from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest


class TestVersion(BaseAPIPortalTest):
    @patch.dict(os.environ, {"COMMIT_SHA": "test"})
    def test_get(self):
        response = self.app.get("/dp/v1/deployed_version")
        self.assertEqual("test", json.loads(response.data)["Data Portal"])
