import json
import os
import unittest
import urllib.parse
import time
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest

from tests.unit.backend.corpora.api_server.mock_auth import TOKEN_EXPIRES


@unittest.skipIf(
    os.environ["DEPLOYMENT_STAGE"] != "test",
    f"Does not run DEPLOYMENT_STAGE:{os.environ['DEPLOYMENT_STAGE']}",
)
class TestAuth(BaseAuthAPITest):
    def check_user_info(self, userinfo):
        self.assertEqual(userinfo["is_authenticated"], True)
        self.assertEqual(userinfo["id"], "test_user_id")
        self.assertEqual(userinfo["name"], "Fake User")
        self.assertEqual(userinfo["email"], "fake_user@email.com")
        self.assertEqual(userinfo["email_verified"], True)

    def check_set_cookie_is_secure(self, cookie, expected_flags):
        actual_flags = cookie.split("; ")
        for expected_flag in expected_flags:
            self.assertIn(expected_flag, actual_flags)

    def test__auth_flow(self):
        headers = dict(host="localhost")

        with self.subTest("userinfo_not_authenticated"):
            response = self.app.get("/dp/v1/userinfo", headers=headers)
            self.assertEqual(401, response.status_code)
            body = json.loads(response.data)
            self.assertEqual(body["detail"], "No authorization token provided")

        with self.subTest("login"):
            response = self.app.get("/dp/v1/login", headers=headers)
            self.assertEqual(response.status_code, 302)
            location = response.headers["Location"]
            split = urllib.parse.urlsplit(location)
            args = dict(urllib.parse.parse_qsl(split.query))
            self.assertTrue(location.startswith(f"{self.auth_config.api_authorize_url}"))
            self.assertIn("response_type=code", location)
            self.assertEqual(args["client_id"], self.auth_config.client_id)
            self.assertEqual(args["response_type"], "code")
            self.assertIn("/dp/v1/oauth2/callback", args["redirect_uri"])
            # Test session cookie
            self.check_set_cookie_is_secure(response.headers["Set-Cookie"], ["HttpOnly", "Secure", "SameSite=Lax"])

            # follow redirect
            test_url = f"/dp/v1/oauth2/callback?code=fakecode&state={args['state']}"
            response = self.app.get(test_url, headers=dict(host="localhost", Cookie=response.headers["Set-Cookie"]))
            self.assertEqual(response.status_code, 302)
            self.assertEqual(response.headers["Location"], self.auth_config.redirect_to_frontend)
            self.assertIn("Set-Cookie", response.headers)
            # Test cxguser cookie
            self.check_set_cookie_is_secure(
                response.headers["Set-Cookie"], expected_flags=["HttpOnly", "Secure", "SameSite=Strict"]
            )
            cxguser_cookie = response.headers["Set-Cookie"]

            # check userinfo
            response = self.app.get("/dp/v1/userinfo", headers=dict(host="localhost", Cookie=cxguser_cookie))
            self.assertEqual(200, response.status_code)

            body = json.loads(response.data)
            self.check_user_info(body)
            self.assertNotIn("Set-Cookie", response.headers)  # no cookie expected

            # sleep so the token expires, then try userinfo again, verify it refreshed
            time.sleep(TOKEN_EXPIRES + 1)

            # check the userinfo again (token has now expired, make sure it is refreshed)
            response = self.app.get("/dp/v1/userinfo", headers=dict(host="localhost", Cookie=cxguser_cookie))
            self.assertEqual(200, response.status_code)
            body = json.loads(response.data)
            self.check_user_info(body)
            self.assertIn("Set-Cookie", response.headers)
            cxguser_cookie = response.headers["Set-Cookie"]

            # check the userinfo again (make sure the replacement cookie works)
            response = self.app.get("/dp/v1/userinfo", headers=dict(host="localhost", Cookie=cxguser_cookie))
            self.assertEqual(200, response.status_code)
            body = json.loads(response.data)
            self.check_user_info(body)
            self.assertNotIn("Set-Cookie", response.headers)

        with self.subTest("login_redirect"):
            response = self.app.get("/dp/v1/login?redirect=?showCC=1", headers=headers)
            self.assertEqual(response.status_code, 302)
            location = response.headers["Location"]
            split = urllib.parse.urlsplit(location)
            args = dict(urllib.parse.parse_qsl(split.query))

            # follow redirect
            test_url = f"/dp/v1/oauth2/callback?code=fakecode&state={args['state']}"
            response = self.app.get(test_url, headers=dict(host="localhost", Cookie=response.headers["Set-Cookie"]))
            self.assertEqual(response.status_code, 302)
            self.assertEqual(response.headers["Location"], self.auth_config.redirect_to_frontend + "?showCC=1")

        with self.subTest("logout"):
            response = self.app.get("/dp/v1/logout", headers=headers)
            self.assertEqual(response.status_code, 302)
            location = response.headers["Location"]
            split = urllib.parse.urlsplit(location)
            args = dict(urllib.parse.parse_qsl(split.query))
            self.assertTrue(location.startswith(f"{self.auth_config.api_base_url}/v2/logout"))
            self.assertEqual(args["returnTo"], self.auth_config.redirect_to_frontend)
            self.assertTrue("Set-Cookie" in response.headers)
            self.assertTrue(response.headers["Set-Cookie"].startswith("cxguser=;"))

            # check userinfo
            response = self.app.get("/dp/v1/userinfo", headers=dict(host="localhost"))
            self.assertEqual(401, response.status_code)
