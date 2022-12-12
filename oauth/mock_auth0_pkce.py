import time

from flask import Flask, request, redirect, make_response, jsonify, session
from flask_cors import CORS
from urllib.parse import quote

# seconds until the token expires
from jose import jwt


TOKEN_EXPIRES = 5 * 60

# A mocked out oauth server, which serves all the endpoints needed by the oauth type.
class MockOauthApp:
    def __init__(self, port, additional_scope=None, token_duration=TOKEN_EXPIRES):
        self.port = port
        self.additional_scope = additional_scope if additional_scope else []
        self.token_duration = int(token_duration)

        # mock flask app
        self.app = Flask("mock_oauth_app")
        self.app.config.update(
            SECRET_KEY="secret_key",
            SESSION_COOKIE_HTTPONLY=True,
            SESSION_COOKIE_SAMESITE="Lax",
            JSON_SORT_KEYS=True,
        )

        self.app.add_url_rule("/authorize", view_func=self.api_authorize)
        self.app.add_url_rule("/oauth/token", view_func=self.api_oauth_token, methods=["OPTIONS", "POST"])
        self.app.add_url_rule("/v2/logout", view_func=self.api_logout)
        self.app.add_url_rule("/userinfo", view_func=self.api_userinfo)
        self.app.add_url_rule("/.well-known/openid-configuration", view_func=self.api_openid_configuration)
        self.app.add_url_rule("/.well-known/jwks.json", view_func=self.api_jwks)
        CORS(self.app, origins=["https://frontend.corporanet.local:3000"], allow_headers=["content-type", "auth0-client"])

    def api_authorize(self):
        session["code_djh"] = request.args.get("code_challenge")
        callback = request.args.get("redirect_uri")
        state = request.args.get("state")
        scope = request.args.get("scope")
        return redirect(callback + f"?code=fakecode&scope={quote(scope)}&state={quote(state)}")

    def api_userinfo(self):
        return dict(
            name="Fake User",
            sub="test_user_id",
            id="test_user_id",
            email="fake_user@email.com",
            email_verified=True,
        )

    def api_oauth_token(self):
        # headers = [
        #     ("Access-Control-Allow-Origin", "https://frontend.corporanet.local:3000"),
        #     ("Access-Control-Allow-Headers", "content-type"),
        #     ("Access-Control-Allow-Headers", "auth0-client"),
        # ]
        if request.method == "OPTIONS":
            return make_response("", 200)
        elif request.method == "POST":
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
            return make_response(jsonify(r), 200)

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
        self.app = MockOauthApp(443, self.additional_scope, self.token_duration).app


def make_token(token_claims: dict, token_duration: int = 5, additional_scope: list = None) -> str:
    additional_scope = additional_scope if additional_scope else []
    expires_at = time.time() + token_duration
    headers = dict(alg="RS256", kid="fake_kid")
    token_claims.update(exp=expires_at, scope=additional_scope)

    token = jwt.encode(claims=token_claims, key="mysecret", algorithm="HS256", headers=headers)
    return token


mock_auth_server = MockOauthServer()
app = mock_auth_server.app
