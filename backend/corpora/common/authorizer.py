import os
from functools import lru_cache

import requests
from chalice import UnauthorizedError
from jose import jwt
from jose.exceptions import ExpiredSignatureError, JWTError, JWTClaimsError

from .corpora_config import CorporaAuthConfig


def assert_authorized_token(token: str) -> dict:
    """
    Determines if the Access Token is valid and return the decoded token. Userinfo is added to the token if it exists.
    :param token: The token
    :return: The decoded access token and userinfo.
    """
    try:
        unverified_header = jwt.get_unverified_header(token)
    except JWTError:
        raise UnauthorizedError(msg="Unable to parse authentication token.")
    auth_config = CorporaAuthConfig()
    auth0_domain = auth_config.internal_url
    audience = auth_config.audience
    public_keys = get_public_keys(auth0_domain)
    public_key = public_keys.get(unverified_header["kid"])
    if public_key:
        algorithms = ["RS256"]
        options = {}
        # in some test situations ignore verifying the signature and issuer
        if os.environ.get("IS_DOCKER_DEV") or (
            os.environ.get("DEPLOYMENT_STAGE") == "test" and (not public_key.get("n") or not public_key.get("e"))
        ):
            options = {"verify_signature": False, "verify_iss": False, "verify_at_hash": False}
        try:
            if not auth0_domain.endswith("/"):
                auth0_domain += "/"
            payload = jwt.decode(
                token, public_key, algorithms=algorithms, audience=audience, issuer=auth0_domain, options=options
            )
        except ExpiredSignatureError:
            raise
        except JWTClaimsError:
            raise UnauthorizedError(msg="Incorrect claims, please check the audience and issuer.")
        except Exception:
            raise UnauthorizedError(msg="Unable to parse authentication token.")

        return payload

    raise UnauthorizedError(msg="Unable to find appropriate key")


def assert_authorized(headers: dict) -> dict:
    """
    Determines if the Access Token is valid and return the decoded token. Userinfo is added to the token if it exists.
    :param headers: The http headers from the request.
    :return: The decoded access token and userinfo.
    """
    try:
        token = get_token_auth_header(headers)
        return assert_authorized_token(token)
    except ExpiredSignatureError:
        raise UnauthorizedError(msg="Token is expired.")


def get_token_auth_header(headers: dict) -> str:
    """Obtains the Access Token from the Authorization Header """

    auth_header = headers.get("Authorization", None)
    if not auth_header:
        raise UnauthorizedError(msg="Authorization header is expected")

    parts = auth_header.split()

    if parts[0].lower() != "bearer":
        raise UnauthorizedError(msg="Authorization header must start with Bearer")
    elif len(parts) == 1:
        raise UnauthorizedError(msg="Token not found")
    elif len(parts) > 2:
        raise UnauthorizedError(msg="Authorization header must be Bearer token")

    token = parts[1]
    return token


def get_userinfo(token: str) -> dict:
    if token is None:
        userinfo = dict(is_authenticated=False)
        return userinfo

    payload = assert_authorized_token(token)

    userinfo = dict(
        is_authenticated=True,
        id=payload.get("sub"),
        name=payload.get("name"),
        email=payload.get("email"),
        email_verified=payload.get("email_verified"),
    )
    return userinfo


@lru_cache(maxsize=32)
def get_openid_config(openid_provider: str):
    """
    :param openid_provider: the openid provider's domain.
    :return: the openid configuration
    """
    res = requests.get("{op}/.well-known/openid-configuration".format(op=openid_provider))
    res.raise_for_status()
    return res.json()


@lru_cache(maxsize=32)
def get_public_keys(openid_provider: str):
    """
    Fetches the public key from an OIDC Identity provider to verify the JWT.
    :param openid_provider: the openid provider's domain.
    :return: Public Keys
    """
    keys = requests.get(get_openid_config(openid_provider)["jwks_uri"]).json()["keys"]
    return {key["kid"]: key for key in keys}
