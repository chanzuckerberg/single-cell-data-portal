import os
from functools import lru_cache

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util import Retry

from backend.corpora.common.corpora_config import CorporaAuthConfig
from backend.corpora.common.utils.http_exceptions import UnauthorizedError
from backend.corpora.common.utils.jwt import jwt_decode, get_unverified_header

auth0_session_with_retry = requests.Session()
retry_config = Retry(total=3, backoff_factor=1, status_forcelist=CorporaAuthConfig().retry_status_forcelist)
auth0_session_with_retry.mount("https://", HTTPAdapter(max_retries=retry_config))


def assert_authorized_token(token: str, audience: str = None) -> dict:
    """
    Determines if the Access Token is valid and return the decoded token. Userinfo is added to the token if it exists.
    :param token: The token
    :return: The decoded access token and userinfo.
    """
    unverified_header = get_unverified_header(token)
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
        payload = jwt_decode(
            token, public_key, algorithms=algorithms, audience=use_audience, issuer=issuer, options=options
        )
        return payload

    raise UnauthorizedError(detail="Unable to find appropriate key")


def get_userinfo_from_auth0(token: str) -> dict:
    auth_config = CorporaAuthConfig()
    res = auth0_session_with_retry.get(auth_config.api_userinfo_url, headers={"Authorization": f"Bearer {token}"})
    res.raise_for_status()
    return res.json()


@lru_cache(maxsize=32)
def get_openid_config(openid_provider: str):
    """
    :param openid_provider: the openid provider's domain.
    :return: the openid configuration
    """
    res = auth0_session_with_retry.get("{op}/.well-known/openid-configuration".format(op=openid_provider))
    res.raise_for_status()
    return res.json()


@lru_cache(maxsize=32)
def get_public_keys(openid_provider: str):
    """
    Fetches the public key from an OIDC Identity provider to verify the JWT.
    :param openid_provider: the openid provider's domain.
    :return: Public Keys
    """
    keys = auth0_session_with_retry.get(get_openid_config(openid_provider)["jwks_uri"]).json()["keys"]
    return {key["kid"]: key for key in keys}
