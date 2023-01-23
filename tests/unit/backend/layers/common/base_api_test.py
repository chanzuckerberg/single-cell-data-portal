import base64
import json
import time
import unittest
from unittest.mock import patch
from backend.layers.thirdparty.cdn_provider_interface import CDNProviderInterface
from tests.unit.backend.api_server.config import TOKEN_EXPIRES
from tests.unit.backend.layers.common.base_test import BaseTest


class BaseAuthAPITest(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.mock_assert_authorized_token = patch(
            "backend.portal.api.app.v1.authentication.assert_authorized_token",
            side_effect=self._mock_assert_authorized_token,
        )
        self.mock_assert_authorized_token.start()

        self.mock_config = patch(
            "backend.curation.api.v1.curation.collections.common.get_collections_base_url",
            return_value="https://frontend.corporanet.local:3000",
        )
        self.mock_config.start()

    def tearDown(self):
        super().tearDown()
        self.mock_assert_authorized_token.stop()

    def make_owner_header(self):
        return {"Authorization": "Bearer " + "owner", "Content-Type": "application/json"}

    def make_super_curator_header(self):
        return {"Authorization": "Bearer " + "super", "Content-Type": "application/json"}

    def make_not_owner_header(self):
        return {"Authorization": "Bearer " + "not_owner", "Content-Type": "application/json"}

    def _mock_assert_authorized_token(self, token: str, audience: str = None):
        if token == "owner":
            return {"sub": "test_user_id", "email": "fake_user@email.com", "scope": []}
        elif token == "not_owner":
            return {"sub": "someone_else", "email": "fake_user@email.com", "scope": []}
        elif token == "super":
            return {"sub": "super", "email": "fake_user@email.com", "scope": ["write:collections"]}
        else:
            raise Exception()


class BaseAPIPortalTest(BaseAuthAPITest, BaseTest):
    def setUp(self):
        super().setUp()

        # TODO: this can be improved, but the current authorization method requires it
        self.mock = patch("backend.common.corpora_config.CorporaAuthConfig.__getattr__", return_value="mock_audience")
        self.mock.start()

        self.cloudfront_provider = CDNProviderInterface()

        from backend.api_server.app import app

        self.app = app.test_client(use_cookies=False)

        # Mock all the dependencies of the API classes
        self.mock_business_logic = patch("backend.portal.api.providers._business_logic", new=self.business_logic)
        self.mock_business_logic.start()

        self.mock_cloudfront_provider = patch(
            "backend.portal.api.providers._cloudfront_provider", new=self.cloudfront_provider
        )
        self.mock_cloudfront_provider.start()

    def tearDown(self):
        super().tearDown()
        # Disable mocking of business logic and cloudfront provider
        self.mock_business_logic.stop()
        self.mock_cloudfront_provider.stop()

    def get_mock_cxguser_token(self, user="owner"):
        """
        Generated an auth token for testing.
        :param user: the type of use the token will simulate.
        :return:
        """
        return f"Bearer {user}"
