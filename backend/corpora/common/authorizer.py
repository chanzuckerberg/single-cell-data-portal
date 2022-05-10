import os
from functools import lru_cache

import requests
from jose import jwt
from jose.exceptions import ExpiredSignatureError, JWTError, JWTClaimsError

from .corpora_config import CorporaAuthConfig
from backend.corpora.common.utils.http_exceptions import UnauthorizedError


def assert_authorized_token(token: str, audience: str = None) -> dict:
    """
    Determines if the Access Token is valid and return the decoded token. Userinfo is added to the token if it exists.
    :param token: The token
    :return: The decoded access token and userinfo.
    """
    try:
        unverified_header = jwt.get_unverified_header(token)
    except JWTError:
        raise UnauthorizedError(detail="Unable to parse authentication token.")
    auth_config = CorporaAuthConfig()
    auth0_domain = auth_config.internal_url
    # If we're using an id_token (for userinfo), we need a difference audience, which gets passed in.
    # Otherwise use auth_config.api_audience
    use_audience = audience or auth_config.api_audience
    public_keys = get_public_keys(auth0_domain)
    public_key = public_keys.get(unverified_header["kid"])
    if public_key:
        algorithms = ["RS256"]
        options = {}
        issuer = auth_config.issuer

        # in some test situations ignore verifying the signature and issuer
        if os.environ.get("IS_DOCKER_DEV") or (
            os.environ.get("DEPLOYMENT_STAGE") == "test" and (not public_key.get("n") or not public_key.get("e"))
        ):
            options = {"verify_signature": False, "verify_iss": False, "verify_at_hash": False}
        try:
            payload = jwt.decode(
                token, public_key, algorithms=algorithms, audience=use_audience, issuer=issuer, options=options
            )
        except ExpiredSignatureError:
            raise
        except JWTClaimsError:
            raise UnauthorizedError(detail="Incorrect claims, please check the audience and issuer.")
        except Exception:
            raise UnauthorizedError(detail="Unable to parse authentication token.")

        return payload

    raise UnauthorizedError(detail="Unable to find appropriate key")


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
        raise UnauthorizedError(detail="Token is expired.")


def get_token_auth_header(headers: dict) -> str:
    """Obtains the Access Token from the Authorization Header"""

    auth_header = headers.get("Authorization", None)
    if not auth_header:
        raise UnauthorizedError(detail="Authorization header is expected")

    parts = auth_header.split()

    if parts[0].lower() != "bearer":
        raise UnauthorizedError(detail="Authorization header must start with Bearer")
    elif len(parts) == 1:
        raise UnauthorizedError(detail="Token not found")
    elif len(parts) > 2:
        raise UnauthorizedError(detail="Authorization header must be Bearer token")

    token = parts[1]
    return token


def get_userinfo_from_auth0(token: str) -> dict:
    auth_config = CorporaAuthConfig()
    res = requests.get(auth_config.api_userinfo_url, headers={"Authorization": f"Bearer {token}"})
    res.raise_for_status()
    return res.json()


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
