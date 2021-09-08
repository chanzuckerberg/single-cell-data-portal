import os
import typing

from backend.corpora.api_server.app import app
from backend.corpora.common.corpora_config import CorporaAuthConfig

from tests.unit.backend.corpora.api_server.mock_auth import MockOauthServer
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class BaseAPITest(DataPortalTestCase):
    """
    Provide access to the test Corpora API. All test for the Corpora API should inherit this class.
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


class BaseAuthAPITest(BaseAPITest):
    @staticmethod
    def get_mock_server_and_auth_config(additional_scope=None):
        mock_oauth_server = MockOauthServer(additional_scope)
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
        }
        auth_config.set(authconfig)
        return (mock_oauth_server, auth_config)

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        (mock_oauth_server, auth_config) = BaseAuthAPITest.get_mock_server_and_auth_config()
        cls.mock_oauth_server = mock_oauth_server
        cls.auth_config = auth_config

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        cls.mock_oauth_server.terminate()


class BasicAuthAPITestCurator(BaseAPITest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        (mock_oauth_server, auth_config) = BaseAuthAPITest.get_mock_server_and_auth_config("write:collections")
        cls.mock_oauth_server = mock_oauth_server
        cls.auth_config = auth_config

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        cls.mock_oauth_server.terminate()
