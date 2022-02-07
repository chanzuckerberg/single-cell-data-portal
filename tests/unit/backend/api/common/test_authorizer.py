import json
import os
import unittest

import requests

from backend.api.data_portal.common.authorizer import assert_authorized
from backend.api.data_portal.config.app_config import AuthConfig
from backend.api.data_portal.common.utils.aws import AwsSecret
from backend.api.data_portal.common.utils.exceptions import UnauthorizedError
from tests.unit.backend.api.fixtures.test_config import TemporaryTestConfigChange

"""
These test may start failing if the monthly allowance of Auth0 machine to machine access tokens is exhasted.
"""

# TODO: THIS IS NOT A UNIT TEST. It calls out to real auth0 server!
# TODO: Do we actually need to set anything in AWS Secrets or can't we just now use DefaultConfigPropSource?
@unittest.skipIf(
    os.environ["DEPLOYMENT_STAGE"] != "test",
    f"Does not run DEPLOYMENT_STAGE:{os.environ['DEPLOYMENT_STAGE']}",
)
class TestAuthorizer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        # TODO: A clearer name for this secret might be "corpora/backend/live_test/auth0-secret"
        secret_name = "corpora/cicd/test/auth0-secret"  # using the same secret for all non production stages.
        cls.auth0_secret = json.loads(AwsSecret(secret_name).value)
        cls.auth0_secret["audience"] = f"https://api.{os.getenv('DEPLOYMENT_STAGE')}.corpora.cziscience.com"

        cls.test_auth_config = dict(
                api_base_url="https://czi-single-cell.auth0.com",
                audience=cls.auth0_secret["audience"],
                client_id=cls.auth0_secret["client_id"])

    @unittest.skipIf(os.getenv("IS_DOCKER_DEV"), "No access to auth0 in docker")
    def test_positive(self):
        with TemporaryTestConfigChange(AuthConfig(), self.test_auth_config):
            token = self.get_auth_token()
            assert_authorized({"Authorization": f"bearer {token['access_token']}"})

    @unittest.skipIf(os.getenv("IS_DOCKER_DEV"), "No access to auth0 in docker")
    def test_not_bearer(self):
        sample_non_bearer_auth_token = {
            "Authorization": "Basic czZCaGRSa3F0MzpnWDFmQmF0M2JW",
            "Content-Type": "application/x-www-form-urlencoded",
            "Accept": "application/json",
        }
        with self.assertRaises(UnauthorizedError):
            assert_authorized(sample_non_bearer_auth_token)

    @unittest.skipIf(os.getenv("IS_DOCKER_DEV"), "No access to auth0 in docker")
    def test_invalid_token(self):
        token = self.get_auth_token()
        header, msg, hash = token["access_token"].split(".")
        with self.subTest("short hash"):
            with self.assertRaises(UnauthorizedError):
                assert_authorized({"Authorization": f"bearer {token['access_token'][:-1]}"})

        with self.subTest("wrong hash"):
            with self.assertRaises(UnauthorizedError):
                _token = ".".join([header, msg, "0" * len(hash)])
                assert_authorized({"Authorization": f"bearer {_token}"})

        with self.subTest("wrong header"):
            with self.assertRaises(UnauthorizedError):
                _token = ".".join(["0" * len(header), msg, hash])
                assert_authorized({"Authorization": f"bearer {_token}"})

        with self.subTest("short header"):
            with self.assertRaises(UnauthorizedError):
                _token = token["access_token"][1:]
                assert_authorized({"Authorization": f"bearer {_token}"})

        with self.subTest("wrong message"):
            with self.assertRaises(UnauthorizedError):
                _token = ".".join([header, "0" * len(msg), hash])
                assert_authorized({"Authorization": f"bearer {_token}"})

    def get_auth_token(self) -> dict:
        return requests.post(
            "https://czi-single-cell.auth0.com/oauth/token",
            json=self.auth0_secret,
            headers={"content-type": "application/json"},
        ).json()
