import base64
import contextlib
import json
import logging
import os
from typing import Optional
from urllib.parse import urlencode

import requests  # type: ignore
from authlib.integrations.flask_client import OAuth
from authlib.integrations.flask_client.remote_app import FlaskRemoteApp
from flask import Response, after_this_request, current_app, g, jsonify, make_response, redirect, request, session

from backend.common.authorizer import assert_authorized_token, get_userinfo_from_auth0
from backend.common.corpora_config import CorporaAuthConfig

# global oauth client
from backend.common.utils.http_exceptions import ExpiredCredentialsError, UnauthorizedError

oauth_client = None


def get_oauth_client(config: CorporaAuthConfig) -> FlaskRemoteApp:
    """Create an oauth client on the first invocation, then return oauth client for subsequent calls.

    :param config:  An object containing the auth configuration.
    :return: The oauth client.
    """

    logging.error(f"DEPLOYMENT_STAGE {os.environ['DEPLOYMENT_STAGE']}")

    global oauth_client
    if oauth_client and os.environ["DEPLOYMENT_STAGE"] != "test":
        # tests may have different configs
        return oauth_client

    code_challenge_method = None
    with contextlib.suppress(RuntimeError):
        code_challenge_method = config.code_challenge_method

    oauth = OAuth(current_app)
    oauth_client = oauth.register(
        "oauth",
        code_challenge_method=code_challenge_method,
        client_id=config.client_id,
        client_secret=config.client_secret,
        api_base_url=config.api_base_url,
        refresh_token_url=config.api_token_url,
        access_token_url=config.api_token_url,
        authorize_url=config.api_authorize_url,
        client_kwargs={"scope": "openid profile email offline_access write:collections"},
        authorize_params={"audience": config.api_audience},
    )
    return oauth_client


def login() -> Response:
    """API call: initiate the login process."""
    config = CorporaAuthConfig()
    redirect = request.args.get("redirect", "")
    if os.getenv("FRONTEND_URL"):
        return_to = f"{os.getenv('FRONTEND_URL')}{redirect}"
    else:
        return_to = f"{config.redirect_to_frontend}{redirect}"
    # save the return path in the session cookie, accessed in the callback function
    session["oauth_corpora_callback_redirect"] = return_to
    client = get_oauth_client(config)

    if os.getenv("API_URL"):
        callbackurl = f"{os.getenv('API_URL')}/dp/v1/oauth2/callback"
    else:
        callbackurl = f"{config.callback_base_url}/dp/v1/oauth2/callback"
    response = client.authorize_redirect(redirect_uri=callbackurl)
    return response


def logout() -> Response:
    """API call: logout of the system."""
    config = CorporaAuthConfig()
    client = get_oauth_client(config)
    params = {"returnTo": config.redirect_to_frontend, "client_id": config.client_id}
    response = redirect(client.api_base_url + "/v2/logout?" + urlencode(params))
    # remove the cookie
    remove_token(config.cookie_name)
    return response  # type: ignore


def oauth2_callback() -> Response:
    """API call: redirect from the auth server after login successful."""
    config = CorporaAuthConfig()
    client = get_oauth_client(config)
    try:
        token = client.authorize_access_token()
        # write the cookie
        save_token(config.cookie_name, token)
    except Exception as e:
        current_app.logger.exception("Unable to authorize access token.")
        # remove the token
        remove_token(config.cookie_name)
        raise UnauthorizedError(detail=f"Response from oauth server not valid. {e}") from None

    return_to = session.pop("oauth_corpora_callback_redirect", "/")
    return redirect(return_to)  # type: ignore


def save_token(cookie_name: str, token: dict) -> None:
    """Save the token, both in the g scope and in a cookie.

    :param cookie_name: The name of the cookie that is updated or created.
    :param token: A dict containing the token information.
    """
    g.token = token

    @after_this_request
    def _save_cookie(response):
        value = base64.b64encode(json.dumps(dict(token)).encode("utf-8"))
        require_secure_cookies = True
        if os.getenv("DEV_MODE_COOKIES"):
            require_secure_cookies = False
        response.set_cookie(
            cookie_name,
            value,
            httponly=True,
            secure=require_secure_cookies,
            samesite="Strict",
            max_age=24 * 60 * 60 * 90,
        )
        return response


def remove_token(cookie_name: str) -> None:
    """Remove the token, both from the g scope and the cookie

    :param cookie_name:  The name of the cookie to be removed.
    """
    g.pop("token", None)

    @after_this_request
    def _remove_cookie(response):
        response.set_cookie(cookie_name, "", expires=0)
        return response


def decode_token(tokenstr: str) -> dict:
    """Return a token dictionary from the string representation of the token.

    :param tokenstr:  The string representation of the token.
    :return: The token dictionary.
    """
    value = base64.b64decode(tokenstr)
    token = json.loads(value)
    return token


def get_token(cookie_name: str) -> dict:
    """Return the token.

    First look in the g scope, then in the cookie.  It is important to put the
    cookie in the g scope to properly handle refresh tokens, where the new token has not yet been written
    into a cookie before it needs to be accessed.

    :param cookie_name:  The name of the cookie that stores the token.
    :return: The token dictionary.
    """
    if "token" in g:
        return g.token

    tokenstr = request.cookies.get(cookie_name)
    g.token = decode_token(tokenstr)  # type: ignore
    return g.token


def check_token(token: dict) -> dict:
    """Check the validity of the token.

    If the token has expired, attempt to refresh the token.
    if valid, return the payload of the id_token, otherwise throw an UnauthorizedError.

    :param token: a dictionary that contains the token information.
    """
    try:
        payload = assert_authorized_token(token.get("access_token"))  # type: ignore
    except ExpiredCredentialsError:
        # attempt to refresh the token
        auth_config = CorporaAuthConfig()
        try:
            token = refresh_expired_token(token)  # type: ignore
            if token is None:
                raise
            payload = assert_authorized_token(token.get("access_token"))  # type: ignore
            # update the cookie with then refreshed token
            save_token(auth_config.cookie_name, token)
        except ExpiredCredentialsError:
            remove_token(auth_config.cookie_name)
            raise

    return payload


def curation_access_token_func(token):
    return assert_authorized_token(token, CorporaAuthConfig().curation_audience)


def curation_access_token_func_lenient(token: str) -> dict:
    """
    Lenient version that allows endpoints to work even if authentication fails.
    Use this for endpoints that also require public access, so if users end up with a bad token,
    they won't be locked out.
    """
    try:
        return assert_authorized_token(token, CorporaAuthConfig().curation_audience)
    except Exception:
        return {}


def apikey_info_func(tokenstr: str, required_scopes: list) -> dict:
    """Function used by connexion in the securitySchemes.

    The return dictionary must contains a "sub" key.

    :params tokenstr:  A string representation of the token
    :params required_scopes: List of required scopes (currently not used).
    :return: The token dictionary.
    """
    token = decode_token(tokenstr)
    payload = check_token(token)
    return payload


def apikey_info_func_lenient(tokenstr: str, required_scopes: list) -> dict:
    """
    Lenient version that allows endpoints to work even if authentication fails.
    Use this for endpoints that also require public access, so if users end up with a bad token,
    they won't be locked out.
    """
    try:
        return apikey_info_func(tokenstr, required_scopes)
    except Exception:
        return {}


def get_userinfo():
    config = CorporaAuthConfig()
    token = get_token(config.cookie_name)
    return get_userinfo_from_auth0(token.get("access_token"))


def userinfo() -> Response:
    """API call: retrieve the user info from the id token stored in the cookie"""
    userinfo = get_userinfo()
    return make_response(jsonify(userinfo))


def refresh_expired_token(token: dict) -> Optional[dict]:
    """Attempt to refresh the expired token.

    If successful, save and return the new token, otherwise return None.

    :param token: a dictionary that contains the token information.
    :return: Either the new token dict, or None.
    """

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
    request = requests.post(auth_config.api_token_url, urlencode(params), headers=headers)
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
