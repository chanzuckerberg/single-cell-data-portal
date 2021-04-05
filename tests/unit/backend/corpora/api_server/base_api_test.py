import os

from backend.corpora.api_server.app import app
from backend.corpora.common.corpora_config import CorporaAuthConfig

from tests.unit.backend.corpora.api_server.mock_auth import MockOauthServer
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class BaseAPITest(DataPortalTestCase):
    """
    Provide access to the test Corpora API. All test for the Corpora API should inherit this class.
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            cls.app = app.test_client()
            cls.maxDiff = None  # Easier to compare json responses.

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()

    @staticmethod
    def remove_timestamps(body: dict) -> dict:
        """
        A helper function to remove timestamps from the response body.
        :param body: The decoded json response body
        :return: The decode json response body with timestamps removed.
        """

        def _remove_timestamps(jrb):
            if not isinstance(jrb, dict):
                return
            jrb.pop("created_at", None)
            jrb.pop("updated_at", None)
            for value in jrb.values():
                if isinstance(value, dict):
                    _remove_timestamps(value)
                elif isinstance(value, list):
                    for list_value in value:
                        _remove_timestamps(list_value)
            return jrb

        return _remove_timestamps(body)


class BaseAuthAPITest(BaseAPITest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.mock_oauth_server = MockOauthServer()
        cls.mock_oauth_server.start()
        assert cls.mock_oauth_server.server_okay

        # Use the CorporaAuthConfig used by the app
        cls.auth_config = CorporaAuthConfig()

        os.environ["API_BASE_URL"] = f"http://localhost:{cls.mock_oauth_server.port}"
        # Overwrite the environment's auth config with our oidc server's config.
        authconfig = {
            "api_base_url": f"http://localhost:{cls.mock_oauth_server.port}",
            "callback_base_url": cls.auth_config.callback_base_url,
            "redirect_to_frontend": cls.auth_config.redirect_to_frontend,
            "client_id": cls.auth_config.client_id,
            "client_secret": cls.auth_config.client_secret,
            "audience": cls.auth_config.audience,
            "cookie_name": cls.auth_config.cookie_name,
        }
        cls.auth_config.set(authconfig)

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        cls.mock_oauth_server.terminate()
