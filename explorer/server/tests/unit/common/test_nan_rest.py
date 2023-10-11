import math
from http import HTTPStatus
from urllib.parse import quote

import server.tests.decode_fbs as decode_fbs
from server.common.config.app_config import AppConfig
from server.tests import FIXTURES_ROOT
from server.tests.unit import BaseTest

VERSION = "v0.2"
BAD_FILTER = {"filter": {"obs": {"annotation_value": [{"name": "xyz"}]}}}


class WithNaNs(BaseTest):
    """Test Case for endpoints"""

    @classmethod
    def setUpClass(cls):
        app_config = AppConfig()
        app_config.update_server_config(multi_dataset__dataroot=FIXTURES_ROOT, app__flask_secret_key="secret")
        super().setUpClass(app_config)
        cls.app.testing = True
        cls.client = cls.app.test_client()

    def setUp(self):
        self.session = self.client
        response = self.session.get("/d/nan.cxg/api/v0.3/s3_uri")
        s3_uri = quote(quote(response.json, safe=""), safe="")
        self.TEST_URL_BASE = f"/s3_uri/{s3_uri}/api/v0.3/"

    def test_initialize(self):
        endpoint = "schema"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        result = self.session.get(url)
        self.assertEqual(result.status_code, HTTPStatus.OK)

    def test_data(self):
        endpoint = "data/var"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        filter = {"filter": {"var": {"index": [[0, 20]]}}}
        header = {"Accept": "application/octet-stream"}
        result = self.session.put(url, headers=header, json=filter)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertTrue(math.isnan(df["columns"][3][3]))

    def test_annotation_obs(self):
        endpoint = "annotations/obs"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertTrue(math.isnan(df["columns"][2][0]))

    def test_annotation_var(self):
        endpoint = "annotations/var"
        url = f"{self.TEST_URL_BASE}{endpoint}"
        header = {"Accept": "application/octet-stream"}
        result = self.session.get(url, headers=header)
        self.assertEqual(result.status_code, HTTPStatus.OK)
        self.assertEqual(result.headers["Content-Type"], "application/octet-stream")
        df = decode_fbs.decode_matrix_FBS(result.data)
        self.assertTrue(math.isnan(df["columns"][2][0]))
