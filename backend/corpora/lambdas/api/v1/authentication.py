from flask import make_response, jsonify, current_app, request, redirect
import json
from authlib.integrations.flask_client import OAuth
from urllib.parse import urlencode
from functools import lru_cache
import base64
from ....common.authorizer import get_userinfo, assert_authorized_token
from ....common.corpora_config import CorporaAuthConfig


@lru_cache(maxsize=1)
def create_oauth_client(config):
    oauth = OAuth(current_app)
    api_base_url = config.api_base_url
    client = oauth.register(
        "oauth",
        client_id=config.client_id,
        client_secret=config.client_secret,
        api_base_url=api_base_url,
        access_token_url=f"{api_base_url}/oauth/token",
        authorize_url=f"{api_base_url}/authorize",
        client_kwargs={"scope": "openid profile email offline"},
    )
    return client


def login():
    """api call,  initiate the login process"""
    config = CorporaAuthConfig()
    client = create_oauth_client(config)
    callbackurl = f"{config.callback_base_url}/v1/oauth2/callback"
    response = client.authorize_redirect(redirect_uri=callbackurl)
    return response


def logout():
    """api call,  logout of the system"""
    config = CorporaAuthConfig()
    client = create_oauth_client(config)
    params = {"returnTo": config.callback_base_url, "client_id": config.client_id}
    response = redirect(client.api_base_url + "/v2/logout?" + urlencode(params))
    # remove the cookie
    response.set_cookie(config.cookie_name, "", expires=0)
    return response


def oauth2_callback():
    """api call,  redirect from the auth server after login successful"""
    config = CorporaAuthConfig()
    client = create_oauth_client(config)
    token = client.authorize_access_token()
    response = redirect(config.callback_base_url)
    # write the cookie
    value = encode_token(dict(token))
    response.set_cookie(config.cookie_name, value, httponly=True, max_age=24 * 60 * 60 * 90)
    return response


def encode_token(token):
    return base64.b64encode(json.dumps(token).encode("utf-8"))


def decode_token(value):
    value = base64.b64decode(value)
    return json.loads(value)


def apikey_info_func(tokenstr, required_scopes):
    """Function used by connexion in the securitySchemes.
    The return dictionary must contains a "sub" key"""
    token = decode_token(tokenstr)
    payload = assert_authorized_token(token.get("id_token"))
    return payload


def userinfo():
    """api call,  retrieve the user info from the id token stored in the cookie"""
    config = CorporaAuthConfig()
    tokenstr = request.cookies.get(config.cookie_name)
    token = decode_token(tokenstr)
    userinfo = get_userinfo(token.get("id_token"))
    return make_response(jsonify(userinfo))
