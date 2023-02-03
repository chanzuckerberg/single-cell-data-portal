import base64
import json
import time
from unittest.mock import DEFAULT, patch

from jose import jws

from backend.common.auth0_manager import auth0_management_session
from tests.unit.backend.layers.api.config import TOKEN_EXPIRES
from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest


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


class TestKeys(BaseAPIPortalTest):
    @classmethod
    def setUpClass(cls):
        auth0_management_session.domain = "localhost"  # setting this so CorporaAuthConfig isn't used in the test
        super().setUpClass()
        cls.api_key_id = "auth0_generate_api_key_id"
        cls.user_name = "test_user_id"
        cls.email = "fake_user@email.com"
        cls.api_key_secret = "a secret_value"

    @patch.multiple(
        "backend.portal.api.app.v1.auth.keys.auth0_management_session",
        get_user_api_key_identity=DEFAULT,
        store_api_key=DEFAULT,
        link_api_key=DEFAULT,
        delete_api_key=DEFAULT,
    )
    @patch.multiple("backend.portal.api.app.v1.auth.keys", CorporaAuthConfig=DEFAULT, get_userinfo=DEFAULT)
    def test__create_key__201(
        self, CorporaAuthConfig, get_userinfo, get_user_api_key_identity, store_api_key, link_api_key, delete_api_key
    ):
        get_user_api_key_identity.return_value = None
        store_api_key.return_value = self.api_key_id
        link_api_key.return_value = None
        delete_api_key.return_value = None
        CorporaAuthConfig().api_key_secret = self.api_key_secret
        get_userinfo.return_value = {"email": "fake_user@email.com"}
        CorporaAuthConfig().days_to_live = 1
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.post("/dp/v1/auth/key", headers=headers)

        self.assertEqual(201, response.status_code)
        body = json.loads(response.data)
        key = body["key"]
        payload = jws.verify(key, self.api_key_secret, algorithms=["HS256"])
        payload = json.loads(payload.decode("utf-8"))
        self.assertEqual(self.user_name, payload["sub"])
        delete_api_key.assert_not_called()
        store_api_key.assert_called_once_with(key, self.email)
        link_api_key.assert_called_once_with(self.user_name, self.api_key_id)
        get_user_api_key_identity.assert_called_once_with(self.user_name)

    def test__create_key__401(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.post("/dp/v1/auth/key", headers=headers)
        self.assertEqual(response.status_code, 401)

    @patch.multiple(
        "backend.portal.api.app.v1.auth.keys.auth0_management_session",
        get_user_api_key_identity=DEFAULT,
        store_api_key=DEFAULT,
        link_api_key=DEFAULT,
        delete_api_key=DEFAULT,
    )
    @patch.multiple("backend.portal.api.app.v1.auth.keys", CorporaAuthConfig=DEFAULT, get_userinfo=DEFAULT)
    def test__regenerate_key__201(
        self, CorporaAuthConfig, get_userinfo, get_user_api_key_identity, store_api_key, link_api_key, delete_api_key
    ):
        CorporaAuthConfig().api_key_secret = self.api_key_secret
        CorporaAuthConfig().days_to_live = 1
        get_userinfo.return_value = {"email": "fake_user@email.com"}
        get_user_api_key_identity.return_value = {"user_id": self.api_key_id}
        delete_api_key.return_value = None
        store_api_key.return_value = self.api_key_id
        link_api_key.return_value = None
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.post("/dp/v1/auth/key", headers=headers)

        self.assertEqual(201, response.status_code)
        body = json.loads(response.data)
        key = body["key"]
        payload = jws.verify(key, self.api_key_secret, algorithms=["HS256"])
        payload = json.loads(payload.decode("utf-8"))
        self.assertEqual(self.user_name, payload["sub"])
        get_user_api_key_identity.assert_called_once_with(self.user_name)
        delete_api_key.assert_called_once()
        store_api_key.assert_called_once_with(key, self.email)
        link_api_key.assert_called_once_with(self.user_name, self.api_key_id)

    @patch.multiple(
        "backend.portal.api.app.v1.auth.keys.auth0_management_session",
        get_user_api_key_identity=DEFAULT,
        delete_api_key=DEFAULT,
    )
    def test__delete_key__202(self, get_user_api_key_identity, delete_api_key):
        get_user_api_key_identity.return_value = {"user_id": self.api_key_id, "username": "ABCDEF"}
        delete_api_key.return_value = None
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete("/dp/v1/auth/key", headers=headers)
        self.assertEqual(response.status_code, 202)

    def test__delete_key__401(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.delete("/dp/v1/auth/key", headers=headers)
        self.assertEqual(response.status_code, 401)

    @patch("backend.portal.api.app.v1.auth.keys.auth0_management_session.get_user_api_key_identity")
    def test__delete_key__404(self, get_user_api_key_identity):
        get_user_api_key_identity.return_value = None
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete("/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)

    @patch("backend.portal.api.app.v1.auth.keys.auth0_management_session.get_user_api_key_identity")
    def test__get_key__200(self, get_user_api_key_identity):
        get_user_api_key_identity.return_value = {"user_id": self.api_key_id}
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.get("/dp/v1/auth/key", headers=headers)
        self.assertEqual(200, response.status_code)
        get_user_api_key_identity.assert_called_once()

    def test__get_key__401(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get("/dp/v1/auth/key", headers=headers)
        self.assertEqual(response.status_code, 401)

    @patch("backend.portal.api.app.v1.auth.keys.auth0_management_session.get_user_api_key_identity")
    def test__get_key__404(self, get_user_api_key_identity):
        get_user_api_key_identity.return_value = None
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.get("/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)
        get_user_api_key_identity.assert_called_once()
