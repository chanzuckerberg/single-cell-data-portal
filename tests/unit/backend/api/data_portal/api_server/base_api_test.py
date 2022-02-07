import typing

from backend.api.app import app
from backend.api.data_portal.config.app_config import AuthConfig
from tests.unit.backend.api.data_portal.api_server.mock_auth import MockOauthServer
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class BaseAPITest(DataPortalTestCase):
    """
    Provide access to the test Corpora API. All test for the Corpora API should inherit this class.
    """

    maxDiff = None  # Easier to compare json responses.

    def setUp(self):
        super().setUp()
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
    """
    Abstract class for test classes that exercise endpoints requiring user authorization.
    This class will start up a mock OAuth server in a separate process. This allows the endpoints to run through their
    authorization process prior to invoking their core logic.
    """

    @classmethod
    def get_mock_oauth_server(cls):
        return MockOauthServer(additional_scope=None, token_duration=0)

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.mock_oauth_server = cls.get_mock_oauth_server()

        # If debugging tests in IDE, you may need to start the mock oauth server from shell, if it doesn't start.
        # Run this command, and swap out mock server start code by (un)commenting the following lines of code.
        # python ./tests/unit/backend/api/data_portal/api_server/mock_auth.py 18000
        # cls.mock_oauth_server.port = 18000
        cls.mock_oauth_server.start()
        assert cls.mock_oauth_server.server_okay

        cls.orig_api_base_url = AuthConfig().api_base_url
        mock_oauth_server_url = f"http://localhost:{cls.mock_oauth_server.port}"
        AuthConfig().api_base_url = mock_oauth_server_url

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        cls.mock_oauth_server.terminate()
        AuthConfig().api_base_url = cls.orig_api_base_url


class BasicAuthAPITestCurator(BaseAuthAPITest):
    @classmethod
    def get_mock_oauth_server(cls):
        return MockOauthServer(additional_scope="write:collections", token_duration=60)
