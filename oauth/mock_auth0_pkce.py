import time

from flask import Flask, request, redirect, make_response, jsonify, session
from flask_cors import CORS
from urllib.parse import quote

# seconds until the token expires
from jose import jwt, jws


TOKEN_EXPIRES = 5 * 60  # seconds (5 minutes)

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
        self.nonce = None  # used by Auth0 for PKCE

        self.app.add_url_rule("/authorize", view_func=self.api_authorize)
        self.app.add_url_rule("/oauth/token", view_func=self.api_oauth_token, methods=["OPTIONS", "POST"])
        self.app.add_url_rule("/v2/logout", view_func=self.api_logout)
        self.app.add_url_rule("/userinfo", view_func=self.api_userinfo)
        self.app.add_url_rule("/.well-known/openid-configuration", view_func=self.api_openid_configuration)
        self.app.add_url_rule("/.well-known/jwks.json", view_func=self.api_jwks)
        CORS(self.app, origins=["https://frontend.corporanet.local:3000"], allow_headers=["content-type", "auth0-client"])

    def api_authorize(self):
        """
        Skip identity provider interaction and go straight to code exchange
        """
        self.nonce = request.args.get("nonce")
        callback = request.args.get("redirect_uri")
        state = request.args.get("state")
        scope = request.args.get("scope")
        code = "fakecode"
        return redirect(callback + f"?code={code}&scope={quote(scope)}&state={quote(state)}")

    def api_userinfo(self):
        return dict(
            name="Fake User",
            sub="test_user_id",
            id="test_user_id",
            email="fake_user@email.com",
            email_verified=True,
        )

    def api_oauth_token(self):
        if request.method == "OPTIONS":
            return make_response("", 200)
        elif request.method == "POST":
            token_claims = dict(
                name="Fake User",
                sub="test_user_id",
                email="fake_user@email.com",
                email_verified=True,
                nonce=self.nonce,
            )
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
        data = dict(jwks_uri=f"https://oidc.corporanet.local:{self.port}/.well-known/jwks.json")
        return make_response(jsonify(data))

    def api_jwks(self):
        data = {
            "alg": "RS256",
            "kty": "RSA",
            "n": "2tdYgAV6t32HTM_uWdHbyE2CstNBjHVp_tU2qh5MQaf5HQk-ag9aVcQjM39odvZK-XBHbWeTWg7V1-dEqmDBZrRinwzRicB41a"
                 "58Q1fyG4_11KXVN0KQJEqEw3PQszV8DD-oQjA7CkdX5M_YIL9ydl-ij03qHdLM5IIv4VTCbszYqPFdeSkwPWiyGU2TmPdxmkir"
                 "YiU2EZ38-8XNKZCbm4B43km08dwyevHYyEGIppyCV4-33Bu38sI7jzuBpCmB3Xa3np8MQhWyg_dVVtO738K0-En2B43nB7v0Rs"
                 "lqlguYvZ_lNS37UXIJNQej2X3VBzi_S4upYE4sDlyGLRGygQ",
            "e": "AQAB"
        }
        return make_response(jsonify(dict(keys=[data])))


class MockOauthServer:
    def __init__(self, additional_scope=None, token_duration=TOKEN_EXPIRES):
        self.process = None
        self.port = None
        self.server_okay = False
        self.additional_scope = additional_scope
        self.token_duration = token_duration
        self.app = MockOauthApp(443, self.additional_scope, self.token_duration).app


def make_token(token_claims: dict, token_duration: int = 5, additional_scope: list = None) -> str:
    additional_scope = additional_scope if additional_scope else ["openid profile email"]
    now = time.time()
    issued_at = now
    expires_at = now + token_duration
    headers = dict(alg="RS256")
    token_claims.update(
        exp=expires_at,
        iat=issued_at,
        # For some reason, Auth0 requires the client id here for the audience value even though the Auth0Provider takes
        # distinct prop values 'audience' and 'clientId'. Without this configuration, the Auth0 React package rejects
        # the JWT.
        aud="local-client-id",
        scope=additional_scope,
        iss="https://oidc.corporanet.local/",
    )

    private_key_jwk = {
        "alg": "RS256",
        "kty": "RSA",
        "n": "2tdYgAV6t32HTM_uWdHbyE2CstNBjHVp_tU2qh5MQaf5HQk-ag9aVcQjM39odvZK-XBHbWeTWg7V1-dEqmDBZrRinwzRicB41a58Q1"
             "fyG4_11KXVN0KQJEqEw3PQszV8DD-oQjA7CkdX5M_YIL9ydl-ij03qHdLM5IIv4VTCbszYqPFdeSkwPWiyGU2TmPdxmkirYiU2EZ38"
             "-8XNKZCbm4B43km08dwyevHYyEGIppyCV4-33Bu38sI7jzuBpCmB3Xa3np8MQhWyg_dVVtO738K0-En2B43nB7v0RslqlguYvZ_lNS"
             "37UXIJNQej2X3VBzi_S4upYE4sDlyGLRGygQ",
        "e": "AQAB",
        "d": "2jCnHk1YQyY3BhCytn8UQKt3SlBzJFXUrq1qaUb4BOYy7A5RWnGgQa7i4e9_-kwqCHU34g7IzZvI_hCpV65MZdgoFCg1qsBqObJUVt"
             "iSnYR1N-V3pjcJfAWIRU9tn6AN5DB71DI-S0tCiPHprQz0VK2ZaIPojn-kpZhfoKxfhxOCGVpbtW3TGdc0JyKWy00U9QN0JP3LFR27"
             "lT4bdm2KSbbkA-kJW_PoHbBlimOf4Zn1LvYlSvap0xd4imJLrdfQV8KYGYKVjRolF811D60HXiATuWQu263qXL0ab66okiA0KkmTiz"
             "CB_-HdFVTy0oJzIV-9rd_-agWeE8MhNiZQ-Q",
        "p": "9MWKffZT0fqspS0LU6_sprhIf7PFbwvvfrtPajx3JUsY-PbOjvK11D0rKEbVFUYmdfHi5SxZHY4xIyNXM4UkDJETycWTdwuB2HP00o"
             "H_Z6fJy5I4Qdg0cnipwrYWujSexvyHblmnCr4peF3eR4HikTQKSvLb4Qv2QQHyF-cjPE8",
        "q": "5OFK0TAxNVQM3ix-L2akGMzzzyfP3eQaOzKsPwld4XJWCPCG86wGG7nGJNzQivUa2oDyKPoK8GwBb-xY0FmSVd7_GFUi79mTrbx4p3"
             "97j-gM5rKNwfHqICZzchbPOcC-RzGZ-y0Ip6-oNFngYuMbtYWEuqxdhbPWQK_vI95EYC8",
        "dp": "MmIfGcKEimpg8zcMZ9OEkOTJ949Xin4YdCu2MxXzKHgG0ehrqD0JdICKy6WY4uIntItvuuCgD1dfge4aWv6C5xtNyXCj5BM4TQfBU"
              "ztzCTEedorIjbUlRpLTzsKQdw-xxx_f-LT3j1yJSL95q5puupjNrZPD8tT1gXgXDxiCxj0",
        "dq": "BrOuHCf8GOKoKc1FuJ-ZyHwf233_8TBfsEIQlxvwGszdRg-889-ioiczbRrmJFt--Mzzyk1gCqAu_pb-FsO4rDkGYTEE_97wGxM8u"
              "TLbbGBQGGU8D2TFBUH6-Wk_kiJZr_3800UDTt2u4DGNh3e89Pi-5TRs3KBicqx0JSm5-NM",
        "qi": "uo13hZqPx-LDoWTioYtrqm5ul-cnKJwanU2alqmQBPKim19H_dOHxEpVw5ZGspqE9ZnCaIxk3igd2wwsAA1WWnqgnr7D7UJECAv3x"
              "-9elPe4DsyDdk56Q0Snp9ITyCO09cofJW_K0qth-PTVlhqZixb0R6jtYSMqdgC7bi8NHnA"
    }

    token = jws.sign(token_claims, key=private_key_jwk, headers=headers, algorithm="RS256")
    # token = jwt.encode(claims=token_claims, key=private_key_jwk, algorithm="HS256", headers=headers)
    return token


mock_auth_server = MockOauthServer()
app = mock_auth_server.app
