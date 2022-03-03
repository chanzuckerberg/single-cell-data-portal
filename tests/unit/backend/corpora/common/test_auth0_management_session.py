import unittest
from mock import MagicMock

import backend.corpora.common.auth0_management_session
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestAuth0ManagementSession(BaseAuthAPITest):
    def test_session_refresh(self):
        session = backend.corpora.common.auth0_management_session.Auth0ManagementSession("http://localhost:5000")
        session.get_auth0_management_token = MagicMock(return_value="Bearer good")
        session.headers["Authorization"] = "Bearer bad"
        response = session.get(self.auth_config.api_base_url + "/test-refresh")
        response.raise_for_status()
        session.get_auth0_management_token.assert_called_once()


if __name__ == "__main__":
    unittest.main()
