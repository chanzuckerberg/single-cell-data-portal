import json
import os
import sys
import unittest
import urllib.parse
import jose.jwt
import time
from flask import Flask, jsonify, make_response, request, redirect
import random
from multiprocessing import Process
from tests.unit.backend.chalice.api_server import BaseAPITest

# Create a mocked out oauth server, which servers all the endpoints needed by the oauth type.
mock_oauth_app = Flask("mock_oauth_app")

# The port that the mock oauth server will listen on
PORT = random.randint(10000, 12000)

# seconds until the token expires
TOKEN_EXPIRES = 2


@mock_oauth_app.route("/authorize")
def authorize():
    callback = request.args.get("redirect_uri")
    state = request.args.get("state")
    return redirect(callback + f"?code=fakecode&state={state}")


@mock_oauth_app.route("/oauth/token", methods=["POST"])
def token():
    expires_at = time.time()
    headers = dict(alg="RS256", kid="fake_kid")
    payload = dict(name="Fake User", sub="fake_id", email="fake_user@email.com", email_verified=True, exp=expires_at)
    jwt = jose.jwt.encode(claims=payload, key="mysecret", algorithm="HS256", headers=headers)
    r = {
        "access_token": f"access-{time.time()}",
        "id_token": jwt,
        "refresh_token": f"random-{time.time()}",
        "scope": "openid profile email offline",
        "expires_in": TOKEN_EXPIRES,
        "token_type": "Bearer",
        "expires_at": expires_at,
    }
    return make_response(jsonify(r))


@mock_oauth_app.route("/v2/logout")
def logout():
    return_to = request.args.get("returnTo")
    return redirect(return_to)


@mock_oauth_app.route("/.well-known/openid-configuration")
def openid_configuration():
    data = dict(jwks_uri=f"http://localhost:{PORT}/.well-known/jwks.json")
    return make_response(jsonify(data))


@mock_oauth_app.route("/.well-known/jwks.json")
def jwks():
    data = dict(
        alg="RS256",
        kty="RSA",
        use="sig",
        kid="fake_kid",
    )
    return make_response(jsonify(dict(keys=[data])))


# function to launch the mock oauth server
def launch_mock_oauth():
    mock_oauth_app.run(port=PORT)


@unittest.skipIf(
    os.environ["DEPLOYMENT_STAGE"] != "test",
    f"Does not run DEPLOYMENT_STAGE:{os.environ['DEPLOYMENT_STAGE']}",
)
class TestAuth(BaseAPITest, unittest.TestCase):
    def setUp(self):
        self.mock_oauth_process = Process(target=launch_mock_oauth)
        self.mock_oauth_process.start()

    def tearDown(self):
        self.mock_oauth_process.terminate()

    def check_user_info(self, userinfo):
        self.assertEqual(userinfo["is_authenticated"], True)
        self.assertEqual(userinfo["id"], "fake_id")
        self.assertEqual(userinfo["name"], "Fake User")
        self.assertEqual(userinfo["email"], "fake_user@email.com")
        self.assertEqual(userinfo["email_verified"], True)

    def test__auth_flow(self):
        old_path = sys.path.copy()

        def restore_path(p):
            sys.path = p

        sys.path.insert(0, os.path.join(self.corpora_api_dir, "chalicelib"))  # noqa
        self.addCleanup(restore_path, old_path)
        from corpora.common.corpora_config import CorporaAuthConfig

        # Use the CorporaAuthConfig used by the chalice app

        self.auth_config = CorporaAuthConfig()
        self.auth_config._config["api_base_url"] = f"http://localhost:{PORT}"
        self.auth_config._config["callback_base_url"] = "http://localhost:5000"
        self.auth_config.update_defaults()

        headers = dict(host="localhost")

        with self.subTest("userinfo_not_authenticated"):
            response = self.app.get("/dp/v1/userinfo", headers=headers)
            self.assertEqual(401, response.status_code)
            body = json.loads(response.body)
            self.assertEqual(body["detail"], "No authorization token provided")

        with self.subTest("login"):
            response = self.app.get("/dp/v1/login", headers=headers)
            self.assertEqual(response.status_code, 302)
            location = response.headers["Location"]
            split = urllib.parse.urlsplit(location)
            args = dict(urllib.parse.parse_qsl(split.query))
            self.assertTrue(location.startswith(f"{self.auth_config.api_authorize_url}"))
            self.assertTrue("response_type=code" in location)
            self.assertEqual(args["client_id"], self.auth_config.client_id)
            self.assertEqual(args["response_type"], "code")
            self.assertTrue("/dp/v1/oauth2/callback" in args["redirect_uri"])

            # follow redirect
            test_url = f"/dp/v1/oauth2/callback?code=fakecode&state={args['state']}"
            response = self.app.get(test_url, headers=dict(host="localhost", Cookie=response.headers["Set-Cookie"]))
            self.assertEqual(response.status_code, 302)
            self.assertEqual(response.headers["Location"], self.auth_config.redirect_to_frontend)
            self.assertTrue("Set-Cookie" in response.headers)
            cxguser_cookie = response.headers["Set-Cookie"]

            # check userinfo
            response = self.app.get("/dp/v1/userinfo", headers=dict(host="localhost", Cookie=cxguser_cookie))
            self.assertEqual(200, response.status_code)
            body = json.loads(response.body)
            self.check_user_info(body)
            self.assertFalse("Set-Cookie" in response.headers)  # no cookie expected

            # sleep so the token expires, then try userinfo again, verify it refreshed
            time.sleep(TOKEN_EXPIRES + 1)

            # check the userinfo again (token has now expired, make sure it is refreshed)
            response = self.app.get("/dp/v1/userinfo", headers=dict(host="localhost", Cookie=cxguser_cookie))
            self.assertEqual(200, response.status_code)
            body = json.loads(response.body)
            self.check_user_info(body)
            self.assertTrue("Set-Cookie" in response.headers)
            cxguser_cookie = response.headers["Set-Cookie"]

            # check the userinfo again (make sure the replacement cookie works)
            response = self.app.get("/dp/v1/userinfo", headers=dict(host="localhost", Cookie=cxguser_cookie))
            self.assertEqual(200, response.status_code)
            body = json.loads(response.body)
            self.check_user_info(body)
            self.assertFalse("Set-Cookie" in response.headers)

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
