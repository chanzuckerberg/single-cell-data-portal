import json
import unittest

import requests
from dcplib.aws_secret import AwsSecret
from dcp_prototype.backend.browser.code.common.authorizer import assert_authorized
from chalice import UnauthorizedError
import os


@unittest.skipIf(os.getenv("DEPLOYMENT_STAGE","test") != "test", "DEPLOYMENT_STAGE not 'test'")
class TestAuthentication(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.auth0_secret = json.loads(AwsSecret("dcp/backend/browser/test/auth0-secret").value)

    def test_postive(self):
        token = self.get_auth_token()
        assert_authorized({"Authorization": f"Bearer {token['access_token']}"})

    def test_not_bearer(self):
        token = self.get_auth_token()
        with self.assertRaises(UnauthorizedError):
            assert_authorized({"Authorization": f"earer {token['access_token']}"})

    def get_auth_token(self) -> dict:
        return requests.post("https://czi-single-cell.auth0.com/oauth/token",
                             json=self.auth0_secret,
                             headers={"content-type": "application/json"}).json()
