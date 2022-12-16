import base64
import json
import os
import time
import typing
from unittest.mock import patch

from backend.api_server.app import app
from backend.common.corpora_config import CorporaAuthConfig
from tests.unit.backend.api_server.config import TOKEN_EXPIRES
from tests.unit.backend.api_server.mock_auth import MockOauthServer
from tests.unit.backend.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.layers.common.base_test import BaseTest


class BaseAPITest(BaseTest):
    """
    Provide access to the test APIs. All tests for APIs should inherit this class.
    """

    maxDiff = None  # Easier to compare json responses.

    def setUp(self):
        super().setUp()
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)

    @staticmethod
    def remove_timestamps(body: dict, remove: typing.List[str] = None) -> dict:
        """
        A helper function to remove timestamps from the response body.
        :param body: The decoded json response body
        :param remove: Additional attributes to remove.
        :return: The decode json response body with timestamps removed.
        """
        defaults = ["created_at", "updated_at"]
        remove_attributes = remove + defaults if remove else defaults

        def _remove_timestamps(jrb, removing):
            if not isinstance(jrb, dict):
                return
            for rm in removing:
                jrb.pop(rm, None)
            for value in jrb.values():
                if isinstance(value, dict):
                    _remove_timestamps(value, removing)
                elif isinstance(value, list):
                    for list_value in value:
                        _remove_timestamps(list_value, removing)
            return jrb

        return _remove_timestamps(body, remove_attributes)


def mock_assert_authorized_token(token: str, audience: str = None):
    if token == "owner":
        return {"sub": "test_user_id", "email": "fake_user@email.com", "scope": []}
    elif token == "not_owner":
        return {"sub": "someone_else", "email": "fake_user@email.com", "scope": []}
    elif token == "super":
        return {"sub": "super", "email": "fake_user@email.com", "scope": ["write:collections"]}
    else:
        raise Exception()


def get_cxguser_token(user="owner"):
    """
    Generated an auth token for testing.
    :param user: the type of use the token will simulate.
    :return:
    """
    cxguser = base64.b64encode(
        json.dumps(
            {
                "access_token": user,
                "refresh_token": f"random-{time.time()}",
                "scope": "openid profile email offline",
                "expires_in": TOKEN_EXPIRES,
                "token_type": "Bearer",
                "expires_at": TOKEN_EXPIRES,
            }
        ).encode("utf8")
    ).decode("utf8")
    return f"cxguser={cxguser}"


class BaseAuthAPITest(BaseAPITest):
    def setUp(self):
        super().setUp()

        # TODO: this can be improved, but the current authorization method requires it
        self.mock = patch("backend.common.corpora_config.CorporaAuthConfig.__getattr__", return_value="mock_audience")
        self.mock.start()

        self.mock_assert_authorized_token = patch(
            "backend.portal.api.app.v1.authentication.assert_authorized_token",
            side_effect=mock_assert_authorized_token,
        )
        self.mock_assert_authorized_token.start()

    def tearDown(self):
        super().tearDown()
        self.mock_assert_authorized_token.stop()

    def make_owner_header(self):
        return {"Authorization": "Bearer " + "owner", "Content-Type": "application/json"}

    def make_super_curator_header(self):
        return {"Authorization": "Bearer " + "super", "Content-Type": "application/json"}

    def make_not_owner_header(self):
        return {"Authorization": "Bearer " + "not_owner", "Content-Type": "application/json"}


class AuthServerAPITest(BaseAPITest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        (mock_oauth_server, auth_config) = cls.get_mock_server_and_auth_config()
        cls.mock_oauth_server = mock_oauth_server
        cls.auth_config = auth_config

    @staticmethod
    def get_mock_server_and_auth_config(additional_scope=None, token_duration=0):
        mock_oauth_server = MockOauthServer(additional_scope, token_duration)
        mock_oauth_server.start()
        assert mock_oauth_server.server_okay

        # Use the CorporaAuthConfig used by the app
        auth_config = CorporaAuthConfig()

        os.environ["API_BASE_URL"] = f"http://localhost:{mock_oauth_server.port}"
        # Overwrite the environment's auth config with our oidc server's config.
        authconfig = {
            "api_base_url": f"http://localhost:{mock_oauth_server.port}",
            "callback_base_url": auth_config.callback_base_url,
            "redirect_to_frontend": auth_config.redirect_to_frontend,
            "client_id": auth_config.client_id,
            "client_secret": auth_config.client_secret,
            "audience": auth_config.audience,
            "api_audience": auth_config.api_audience,
            "cookie_name": auth_config.cookie_name,
            "auth0_domain": f"localhost:{mock_oauth_server.port}",
            "curation_audience": auth_config.audience,
        }
        auth_config.set(authconfig)
        return (mock_oauth_server, auth_config)

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        cls.mock_oauth_server.terminate()
