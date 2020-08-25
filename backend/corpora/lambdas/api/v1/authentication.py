from flask import make_response, jsonify, current_app, request, redirect
from authlib.integrations.flask_client import OAuth
from urllib.parse import urlencode
from functools import lru_cache
from ....common.authorizer import get_userinfo
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
        client_kwargs={"scope": "openid profile email"},
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
    response.set_cookie(config["cookie_name"], "", expires=0)
    return response


def oauth2_callback():
    """api call,  redirect from the auth server after login successful"""
    config = CorporaAuthConfig()
    client = create_oauth_client(config)
    token = client.authorize_access_token()
    id_token = token.get("id_token")
    response = redirect(config.callback_base_url)
    # write the cookie
    response.set_cookie(config.cookie_name, id_token, httponly=True, max_age=24 * 60 * 60 * 90)
    return response


def userinfo():
    """api call,  retrieve the user info from the id token stored in the cookie"""
    config = CorporaAuthConfig()
    token = request.cookies.get(config.cookie_name)
    userinfo = get_userinfo(token)
    return make_response(jsonify(userinfo))
