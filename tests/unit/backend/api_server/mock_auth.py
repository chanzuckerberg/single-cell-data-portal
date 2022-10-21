import time
import random
import sys

import requests
from flask import Flask, request, redirect, make_response, jsonify
from multiprocessing import Process

# seconds until the token expires
from jose import jwt
from tests.unit.backend.api_server.config import TOKEN_EXPIRES


# A mocked out oauth server, which serves all the endpoints needed by the oauth type.
class MockOauthApp:
    def __init__(self, port, additional_scope=None, token_duration=TOKEN_EXPIRES):
        self.port = port
        self.additional_scope = additional_scope if additional_scope else []
        self.token_duration = int(token_duration)

        # mock flask app
        self.app = Flask("mock_oauth_app")

        self.app.add_url_rule("/authorize", view_func=self.api_authorize)
        self.app.add_url_rule("/oauth/token", view_func=self.api_oauth_token, methods=["POST"])
        self.app.add_url_rule("/v2/logout", view_func=self.api_logout)
        self.app.add_url_rule("/userinfo", view_func=self.api_userinfo)
        self.app.add_url_rule("/.well-known/openid-configuration", view_func=self.api_openid_configuration)
        self.app.add_url_rule("/.well-known/jwks.json", view_func=self.api_jwks)

    def api_authorize(self):
        callback = request.args.get("redirect_uri")
        state = request.args.get("state")
        return redirect(callback + f"?code=fakecode&state={state}")

    def api_userinfo(self):
        return dict(
            name="Fake User",
            sub="test_user_id",
            id="test_user_id",
            email="fake_user@email.com",
            email_verified=True,
        )

    def api_oauth_token(self):
        token_claims = dict(name="Fake User", sub="test_user_id", email="fake_user@email.com", email_verified=True)
        jwt = make_token(token_claims, self.token_duration, self.additional_scope)
        r = {
            "access_token": jwt,
            "id_token": jwt,
            "refresh_token": f"random-{time.time()}",
            "scope": "openid profile email offline " + " ".join(self.additional_scope),
            "expires_in": TOKEN_EXPIRES,
            "token_type": "Bearer",
            "expires_at": self.token_duration,
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
    def __init__(self, additional_scope=None, token_duration=0):
        self.process = None
        self.port = None
        self.server_okay = False
        self.additional_scope = additional_scope
        self.token_duration = token_duration

    def start(self):
        self.port = random.randint(10000, 20000)
        params = [sys.executable, __file__, str(self.port)]
        if self.additional_scope:
            params.append(self.additional_scope)
        if self.token_duration:
            params.append(str(self.token_duration))
        mock_app = MockOauthApp(self.port, self.additional_scope, self.token_duration).app
        self.process = Process(target=mock_app.run, kwargs=dict(port=self.port, debug=True, use_reloader=False))
        self.process.start()
        # Verify that the mock oauth server is ready (accepting requests) before starting the tests.
        self.server_okay = False
        for _ in range(5):
            try:
                response = requests.get(f"http://localhost:{self.port}/.well-known/jwks.json")
                response.raise_for_status()
            except Exception:
                time.sleep(0.1)
            else:
                self.server_okay = True
                break

    def terminate(self):
        self.process.terminate()
        self.process.join()


def make_token(token_claims: dict, token_duration: int = 5, additional_scope: list = None) -> str:
    additional_scope = additional_scope if additional_scope else []
    expires_at = time.time() + token_duration
    headers = dict(alg="RS256", kid="fake_kid")
    token_claims.update(exp=expires_at, scope=additional_scope)

    token = jwt.encode(claims=token_claims, key="mysecret", algorithm="HS256", headers=headers)
    return token
