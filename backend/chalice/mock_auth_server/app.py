import os

import jose.jwt
import time
from chalice import Chalice, Response

app = Chalice(app_name="mock_auth_server")


@app.route("/authorize")
def api_authorize():
    callback = app.current_request.query_params.get("redirect_uri")
    state = app.current_request.query_params.get("state")
    return Response(headers={"Location": f"{callback}?code=fakecode&state={state}"}, status_code=302)


@app.route("/oauth/token")
def api_oauth_token():
    expires_at = time.time()
    headers = dict(alg="RS256", kid="fake_kid")
    payload = dict(
        name="Fake User", sub="test_user_id", email="fake_user@email.com", email_verified=True, exp=expires_at
    )

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
    return r


@app.route("/v2/logout")
def api_logout():
    return_to = app.current_request.query_params.get("returnTo")
    return Response(headers={"Location": return_to}, status_code=302)


@app.route("/.well-known/openid-configuration")
def api_openid_configuration():
    data = dict(jwks_uri=f"{os.getenv('DOMAIN')}/.well-known/jwks.json")
    return data


@app.route("/.well-known/jwks.json")
def api_jwks():
    data = dict(
        alg="RS256",
        kty="RSA",
        use="sig",
        kid="fake_kid",
    )
    return dict(keys=[data])


# seconds until the token expires
TOKEN_EXPIRES = 3600
