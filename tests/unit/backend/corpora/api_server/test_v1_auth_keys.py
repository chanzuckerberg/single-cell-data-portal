import json

from jose import jws
from mock import patch, DEFAULT

from backend.corpora.common.auth0_manager import auth0_management_session
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token


class TestKeys(BaseAuthAPITest):
    @classmethod
    def setUpClass(cls):
        auth0_management_session.domain = "localhost"  # setting this so CorporaAuthConfig isn't used in the test
        super().setUpClass()
        cls.api_key_id = "auth0_generate_api_key_id"
        cls.user_name = "test_user_id"
        cls.email = "fake_user@email.com"
        cls.api_key_secret = "a secret_value"

    @patch.multiple(
        "backend.corpora.lambdas.api.v1.auth.keys.auth0_management_session",
        get_user_api_key_identity=DEFAULT,
        store_api_key=DEFAULT,
        link_api_key=DEFAULT,
        delete_api_key=DEFAULT,
    )
    @patch("backend.corpora.lambdas.api.v1.auth.keys.CorporaAuthConfig")
    def test__create_key__202(
        self, CorporaAuthConfig, get_user_api_key_identity, store_api_key, link_api_key, delete_api_key
    ):
        get_user_api_key_identity.return_value = None
        store_api_key.return_value = self.api_key_id
        link_api_key.return_value = None
        delete_api_key.return_value = None
        CorporaAuthConfig.api_key_secret = self.api_key_secret
        CorporaAuthConfig.days_to_live = 1
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.post("/dp/v1/auth/key", headers=headers)

        self.assertEqual(202, response.status_code)
        body = json.loads(response.data)
        key = body["key"]
        key_name = key.split(".")[-1]
        payload = jws.verify(key, self.api_key_secret, algorithms=["HS256"])
        payload = json.loads(payload.decode("utf-8"))
        self.assertEqual(self.user_name, payload["sub"])
        delete_api_key.assert_not_called()
        store_api_key.assert_called_once_with(key_name, key, self.email)
        link_api_key.assert_called_once_with(self.user_name, self.api_key_id)
        get_user_api_key_identity.assert_called_once_with(self.user_name)

    def test__create_key__401(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.post("/dp/v1/auth/key", headers=headers)
        self.assertEqual(response.status_code, 401)

    @patch.multiple(
        "backend.corpora.lambdas.api.v1.auth.keys.auth0_management_session",
        get_user_api_key_identity=DEFAULT,
        store_api_key=DEFAULT,
        link_api_key=DEFAULT,
        delete_api_key=DEFAULT,
    )
    @patch("backend.corpora.lambdas.api.v1.auth.keys.CorporaAuthConfig")
    def test__regenerate_key__202(
        self, CorporaAuthConfig, get_user_api_key_identity, store_api_key, link_api_key, delete_api_key
    ):
        CorporaAuthConfig.api_key_secret = self.api_key_secret
        CorporaAuthConfig.days_to_live = 1
        get_user_api_key_identity.return_value = {"user_id": self.api_key_id, "username": "abcdefg"}
        delete_api_key.return_value = None
        store_api_key.return_value = self.api_key_id
        link_api_key.return_value = None
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.post("/dp/v1/auth/key", headers=headers)

        self.assertEqual(202, response.status_code)
        body = json.loads(response.data)
        key = body["key"]
        key_name = key.split(".")[-1]
        payload = jws.verify(key, self.api_key_secret, algorithms=["HS256"])
        payload = json.loads(payload.decode("utf-8"))
        self.assertEqual(self.user_name, payload["sub"])
        get_user_api_key_identity.assert_called_once_with(self.user_name)
        delete_api_key.assert_called_once()
        store_api_key.assert_called_once_with(key_name, key, self.email)
        link_api_key.assert_called_once_with(self.user_name, self.api_key_id)

    @patch.multiple(
        "backend.corpora.lambdas.api.v1.auth.keys.auth0_management_session",
        get_user_api_key_identity=DEFAULT,
        delete_api_key=DEFAULT,
    )
    def test__delete_key__201(self, get_user_api_key_identity, delete_api_key):
        get_user_api_key_identity.return_value = {"user_id": self.api_key_id, "username": "ABCDEF"}
        delete_api_key.return_value = None
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete("/dp/v1/auth/key", headers=headers)
        self.assertEqual(response.status_code, 201)

    def test__delete_key__401(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.delete("/dp/v1/auth/key", headers=headers)
        self.assertEqual(response.status_code, 401)

    @patch("backend.corpora.lambdas.api.v1.auth.keys.auth0_management_session.get_user_api_key_identity")
    def test__delete_key__404(self, get_user_api_key_identity):
        get_user_api_key_identity.return_value = None
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete("/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)

    @patch("backend.corpora.lambdas.api.v1.auth.keys.auth0_management_session.get_user_api_key_identity")
    def test__get_key__200(self, get_user_api_key_identity):
        get_user_api_key_identity.return_value = {"user_id": self.api_key_id, "username": "ABCDEF"}
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get("/dp/v1/auth/key", headers=headers)
        print(response.data)
        self.assertEqual(200, response.status_code)
        self.assertEqual(json.loads(response.data), {"id": "ABCDEF"})
        get_user_api_key_identity.assert_called_once()

    def test__get_key__401(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get("/dp/v1/auth/key", headers=headers)
        self.assertEqual(response.status_code, 401)

    @patch("backend.corpora.lambdas.api.v1.auth.keys.auth0_management_session.get_user_api_key_identity")
    def test__get_key__404(self, get_user_api_key_identity):
        get_user_api_key_identity.return_value = None
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get("/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)
        get_user_api_key_identity.assert_called_once()
