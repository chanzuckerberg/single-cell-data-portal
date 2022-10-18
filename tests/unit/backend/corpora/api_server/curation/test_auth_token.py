import unittest
from unittest.mock import patch, Mock

from backend.common.utils.api_key import generate
from tests.unit.backend.corpora.api_server.base_api_test import BaseAPITest


class TestAuthToken(BaseAPITest):
    @patch("backend.portal.api.v1.curation.auth.token.CorporaAuthConfig")
    @patch("backend.portal.api.v1.curation.auth.token.auth0_management_session")
    def test__post_token__201(self, auth0_management_session: Mock, CorporaAuthConfig: Mock):
        test_secret = "password1234"
        test_email = "user@email.com"
        test_user_id = "test_user_id"
        CorporaAuthConfig().api_key_secret = test_secret
        auth0_management_session.get_user_api_key_identity = Mock(return_value={"profileData": {"email": test_email}})
        auth0_management_session.generate_access_token = Mock(return_value={"access_token": "OK"})
        user_api_key = generate(test_user_id, test_secret)
        response = self.app.post("/curation/v1/auth/token", headers={"x-api-key": user_api_key})
        self.assertEqual(201, response.status_code)
        token = response.json["access_token"]
        self.assertEqual("OK", token)
        auth0_management_session.get_user_api_key_identity.assert_called_once_with(test_user_id)

    @patch("backend.portal.api.v1.curation.auth.token.CorporaAuthConfig")
    def test__post_token__401(self, CorporaAuthConfig):
        test_secret = "password1234"
        test_user_id = "test_user_id"
        CorporaAuthConfig().api_key_secret = test_secret
        user_api_key = generate(test_user_id, "not the right secret")
        response = self.app.post("/curation/v1/auth/token", headers={"x-api-key": user_api_key})
        self.assertEqual(401, response.status_code)

    @patch("backend.portal.api.v1.curation.auth.token.CorporaAuthConfig")
    @patch("backend.portal.api.v1.curation.auth.token.auth0_management_session")
    def test__post_token__404(self, auth0_management_session, CorporaAuthConfig):
        test_secret = "password1234"
        test_user_id = "test_user_id"
        CorporaAuthConfig().api_key_secret = test_secret
        auth0_management_session.get_user_api_key_identity = Mock(return_value=None)
        user_api_key = generate(test_user_id, test_secret)
        response = self.app.post("/curation/v1/auth/token", headers={"x-api-key": user_api_key})
        self.assertEqual(404, response.status_code)


if __name__ == "__main__":
    unittest.main()
