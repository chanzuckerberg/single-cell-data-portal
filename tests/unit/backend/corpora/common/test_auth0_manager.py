import unittest
from mock import MagicMock

from backend.corpora.common.auth0_manager import session
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestSession(BaseAuthAPITest):
    def test_session_refresh(self):
        session.domain = "http://localhost:5000"
        session.get_auth0_management_token = MagicMock(return_value="Bearer good")
        session.headers["Authorization"] = "Bearer bad"
        response = session.session.get(self.auth_config.api_base_url + "/test-refresh")
        response.raise_for_status()
        session.get_auth0_management_token.assert_called_once()


if __name__ == "__main__":
    unittest.main()
