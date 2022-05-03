import unittest
from unittest.mock import MagicMock

from backend.corpora.common.auth0_manager import auth0_management_session
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestSession(BaseAuthAPITest):
    def test_session_refresh(self):
        auth0_management_session.domain = "http://localhost:5050"
        auth0_management_session.get_auth0_management_token = MagicMock(return_value="Bearer good")
        auth0_management_session.headers["Authorization"] = "Bearer bad"
        response = auth0_management_session.session.get(self.auth_config.api_base_url + "/test-refresh")
        response.raise_for_status()
        auth0_management_session.get_auth0_management_token.assert_called_once()


if __name__ == "__main__":
    unittest.main()
