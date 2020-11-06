import random
import urllib

import jose.jwt
import time

from flask import Flask, request, redirect, make_response, jsonify

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
    payload = dict(name="Fake User", sub="test_user_id", email="fake_user@email.com", email_verified=True, exp=expires_at)
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


def get_auth_token(app):
    """
    Generated an auth token for testing.
    :param app: a chalice app.
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
