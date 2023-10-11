import json
from http import HTTPStatus

from server.common.config.app_config import AppConfig
from server.tests import FIXTURES_ROOT
from server.tests.unit.common.apis.test_api_v3 import BaseTest

VERSION = "v0.2"


class CorporaRESTAPITest(BaseTest):
    """Confirm endpoints reflect Corpora-specific features"""

    @classmethod
    def setUpClass(cls, app_config=None):
        if not app_config:
            app_config = AppConfig()
        app_config.update_server_config(multi_dataset__dataroot=FIXTURES_ROOT, app__flask_secret_key="test")

        super().setUpClass(app_config)
        cls.app.testing = True
        cls.client = cls.app.test_client()

    def test_config(self):

        test_s3_uri_encoded = self.encode_s3_uri(f"{FIXTURES_ROOT}/schema_2_0_0.cxg")
        url = f"/s3_uri/{test_s3_uri_encoded}/api/v0.3/config"
        header = {"Content-Type": "application/json"}
        result = self.client.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/json")

        result_data = json.loads(result.data)
        self.assertIsInstance(result_data["config"]["corpora_props"], dict)
        self.assertIsInstance(result_data["config"]["parameters"], dict)

        corpora_props = result_data["config"]["corpora_props"]
        result_data["config"]["parameters"]

        self.assertEqual(corpora_props["schema_version"], "2.0.0")

        self.assertEqual(
            corpora_props["title"],
            "Spatiotemporal analysis of human intestinal development at single-cell resolution: Fetal A7",
        )
        # TODO: reinstate? would need to add `default_embedding` to test.cxg file
        # self.assertEqual(parameters["default_embedding"], "tsne")
