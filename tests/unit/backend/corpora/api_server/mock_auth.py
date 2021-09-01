import urllib
import jose.jwt
import time
import random
import sys
import requests
from flask import Flask, request, redirect, make_response, jsonify
import subprocess

# seconds until the token expires
TOKEN_EXPIRES = 2


# A mocked out oauth server, which serves all the endpoints needed by the oauth type.
class MockOauthApp:
    def __init__(self, port, additional_scope = None):
        self.port = port
        self.additional_scope = additional_scope

        # mock flask app
        self.app = Flask("mock_oauth_app")

        self.app.add_url_rule("/authorize", view_func=self.api_authorize)
        self.app.add_url_rule("/oauth/token", view_func=self.api_oauth_token, methods=["POST"])
        self.app.add_url_rule("/v2/logout", view_func=self.api_logout)
        self.app.add_url_rule("/.well-known/openid-configuration", view_func=self.api_openid_configuration)
        self.app.add_url_rule("/.well-known/jwks.json", view_func=self.api_jwks)

    def api_authorize(self):
        callback = request.args.get("redirect_uri")
        state = request.args.get("state")
        return redirect(callback + f"?code=fakecode&state={state}")

    def api_oauth_token(self):
        expires_at = time.time()
        headers = dict(alg="RS256", kid="fake_kid")
        payload = dict(
            name="Fake User", sub="test_user_id", email="fake_user@email.com", email_verified=True, exp=expires_at, scope=self.additional_scope
        )

        jwt = jose.jwt.encode(claims=payload, key="mysecret", algorithm="HS256", headers=headers)
        r = {
            "access_token": jwt,
            "id_token": jwt,
            "refresh_token": f"random-{time.time()}",
            "scope": "openid profile email offline read:collections",
            "expires_in": TOKEN_EXPIRES,
            "token_type": "Bearer",
            "expires_at": expires_at,
        }
        return make_response(jsonify(r))

    def api_logout(self):
        return_to = request.args.get("returnTo")
        return redirect(return_to)

    def api_openid_configuration(self):
        data = dict(jwks_uri=f"http://localhost:{self.port}/.well-known/jwks.json")
        return make_response(jsonify(data))

    def api_jwks(self):
        data = dict(
            alg="RS256",
            kty="RSA",
            use="sig",
            kid="fake_kid",
        )
        return make_response(jsonify(dict(keys=[data])))


class MockOauthServer:
    def __init__(self, additional_scope = None):
        self.process = None
        self.port = None
        self.server_okay = False
        self.additional_scope = additional_scope

    def start(self):
        self.port = random.randint(10000, 20000)
        params = [sys.executable, __file__, str(self.port)]
        if self.additional_scope:
            params.append(self.additional_scope)
        self.process = subprocess.Popen(params)
        # Verify that the mock oauth server is ready (accepting requests) before starting the tests.
        self.server_okay = False
        for _ in range(5):
            try:
                response = requests.get(f"http://localhost:{self.port}/.well-known/jwks.json")
                if response.status_code == 200:
                    self.server_okay = True
                    break
            except Exception:
                pass

            # wait one second and try again
            time.sleep(1)

    def terminate(self):
        self.process.terminate()


def get_auth_token(app):
    """
    Generated an auth token for testing.
    :param app: a WSGI app.
    :return:
    """
    headers = dict(host="localhost")
    response = app.get("/dp/v1/login", headers=headers)
    location = response.headers["Location"]
    split = urllib.parse.urlsplit(location)
    args = dict(urllib.parse.parse_qsl(split.query))

    # follow redirect
    url = f"/dp/v1/oauth2/callback?code=fakecode&state={args['state']}"
    response = app.get(url, headers=dict(host="localhost", Cookie=response.headers["Set-Cookie"]))
    return response.headers["Set-Cookie"]


if __name__ == "__main__":
    port = int(sys.argv[1])
    mock_app = MockOauthApp(port, *sys.argv[2:])
    mock_app.app.run(port=port, debug=True)
