from flask import make_response, jsonify, current_app, request, redirect, after_this_request, g
import json
import os
import requests
from authlib.integrations.flask_client import OAuth
from urllib.parse import urlencode
import base64
from ....common.authorizer import get_userinfo, assert_authorized_token
from backend.corpora.common.corpora_config import CorporaAuthConfig
from jose.exceptions import ExpiredSignatureError
from chalice import UnauthorizedError

# global oauth client
oauth_client = None


def get_oauth_client(config):
    """Create an oauth client on the first invocation, then return oauth client for subsequent calls"""
    global oauth_client
    if oauth_client and os.environ["DEPLOYMENT_STAGE"] != "test":
        # tests may have different configs
        return oauth_client

    oauth = OAuth(current_app)
    api_base_url = config.api_base_url
    oauth_client = oauth.register(
        "oauth",
        client_id=config.client_id,
        client_secret=config.client_secret,
        api_base_url=api_base_url,
        refresh_token_url=f"{api_base_url}/oauth/token",
        access_token_url=f"{api_base_url}/oauth/token",
        authorize_url=f"{api_base_url}/authorize",
        client_kwargs={"scope": "openid profile email offline_access"},
    )
    return oauth_client


def login():
    """api call,  initiate the login process"""
    config = CorporaAuthConfig()
    client = get_oauth_client(config)
    callbackurl = f"{config.callback_base_url}/v1/oauth2/callback"
    response = client.authorize_redirect(redirect_uri=callbackurl)
    return response


def logout():
    """api call,  logout of the system"""
    config = CorporaAuthConfig()
    client = get_oauth_client(config)
    params = {"returnTo": config.callback_base_url, "client_id": config.client_id}
    response = redirect(client.api_base_url + "/v2/logout?" + urlencode(params))
    # remove the cookie
    remove_token(config.cookie_name)
    return response


def oauth2_callback():
    """api call,  redirect from the auth server after login successful"""
    config = CorporaAuthConfig()
    client = get_oauth_client(config)
    token = client.authorize_access_token()
    # write the cookie
    save_token(config.cookie_name, token)
    return redirect(config.callback_base_url)


def save_token(cookie_name, token):
    """Save the token, both in the g scope and in a cookie"""
    g.token = token

    @after_this_request
    def _save_cookie(response):
        value = base64.b64encode(json.dumps(dict(token)).encode("utf-8"))
        response.set_cookie(cookie_name, value, httponly=True, max_age=24 * 60 * 60 * 90)
        return response


def remove_token(cookie_name):
    """Remove the token, both from the g scope and the cookie"""
    g.pop("token", None)

    @after_this_request
    def _remove_cookie(response):
        response.set_cookie(cookie_name, "", expires=0)
        return response


def decode_token(tokenstr):
    """Return a token dictionary from the string representation of the token"""
    value = base64.b64decode(tokenstr)
    token = json.loads(value)
    return token


def get_token(cookie_name):
    """Return the token.  first look in the g scope, then in the cookie.  It is important to put the
    cookie in the g scope to properly handle refresh tokens, where the new token has not yet been written
    into a cookie before it needs to be accessed."""
    if "token" in g:
        return g.token

    tokenstr = request.cookies.get(cookie_name)
    g.token = decode_token(tokenstr)
    return g.token


def check_token(token):
    """check the validity of the token.  If the token has expired, attempt to refresh the token.
    if valid, return the payload of the id_token, otherwise throw an UnauthorizedError"""
    try:
        payload = assert_authorized_token(token.get("id_token"))
        return payload
    except ExpiredSignatureError:
        # attempt to refresh the token
        auth_config = CorporaAuthConfig()
        try:
            token = refresh_expired_token(token)
            if token is None:
                raise
            payload = assert_authorized_token(token.get("id_token"))
            # update the cookie with then refreshed token
            save_token(auth_config.cookie_name, token)
            return payload
        except ExpiredSignatureError:
            remove_token(auth_config.cookie_name)
            raise UnauthorizedError("token is expired")


def apikey_info_func(tokenstr, required_scopes):
    """Function used by connexion in the securitySchemes.
    The return dictionary must contains a "sub" key"""
    token = decode_token(tokenstr)
    payload = check_token(token)
    return payload


def userinfo():
    """api call,  retrieve the user info from the id token stored in the cookie"""
    config = CorporaAuthConfig()
    token = get_token(config.cookie_name)
    userinfo = get_userinfo(token.get("id_token"))
    return make_response(jsonify(userinfo))


def refresh_expired_token(token):
    """Attempt to refresh the expired token.  If successful save return the new token, otherwise return None"""
    auth_config = CorporaAuthConfig()
    refresh_token = token.get("refresh_token")
    if not refresh_token:
        return None
    params = {
        "grant_type": "refresh_token",
        "client_id": auth_config.client_id,
        "refresh_token": refresh_token,
        "client_secret": auth_config.client_secret,
    }
    headers = {"content-type": "application/x-www-form-urlencoded"}
    request = requests.post(f"{auth_config.api_base_url}/oauth/token", urlencode(params), headers=headers)
    if request.status_code != 200:
        # unable to refresh the token
        return None
    data = request.json()
    token = dict(
        access_token=data.get("access_token"),
        id_token=data.get("id_token"),
        refresh_token=data.get("refresh_token", refresh_token),
        expires_at=data.get("expires_at"),
    )
    return token
