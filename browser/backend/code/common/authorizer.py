import os
from functools import lru_cache

import requests
from chalice import UnauthorizedError
from jose import jwt, JWTError

USERINFO_ENDPOINT = "https://czi-single-cell.auth0.com/userinfo"

AUTH0_DOMAIN = "czi-single-cell.auth0.com"
API_AUDIENCE = f"https://api.corpora.{os.environ['DEPLOYMENT_STAGE']}.single-cell.czi.technology"
ALGORITHMS = ["RS256"]


def assert_authorized(headers: dict) -> dict:
    """
    Determines if the Access Token is valid and return the decoded token. Userinfo is
    added to the token if it exists.
    :param headers: The http headers from the request.
    :return: The decoded access token and userinfo.
    """
    token = get_token_auth_header(headers)
    try:
        unverified_header = jwt.get_unverified_header(token)
    except JWTError:
        raise UnauthorizedError(msg="Unable to parse authentication token.")
    public_keys = get_public_keys(AUTH0_DOMAIN)
    public_key = public_keys.get(unverified_header["kid"])
    if public_key:
        try:
            payload = jwt.decode(
                token, public_key, algorithms=ALGORITHMS, audience=API_AUDIENCE, issuer=f"https://{AUTH0_DOMAIN}/"
            )
        except jwt.ExpiredSignatureError:
            raise UnauthorizedError(msg="token is expired")
        except jwt.JWTClaimsError:
            raise UnauthorizedError(msg="incorrect claims, please check the audience and issuer")
        except Exception:
            raise UnauthorizedError(msg="Unable to parse authentication token.")

        if os.getenv("DEPLOYMENT_STAGE", "test").lower() not in ["test", "dev"]:
            payload["userinfo"] = get_userinfo()

        return payload
    raise UnauthorizedError(msg="Unable to find appropriate key")


def get_token_auth_header(headers: dict) -> str:
    """Obtains the Access Token from the Authorization Header
    """
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


def get_userinfo(token):
    resp = requests.post(USERINFO_ENDPOINT, headers={"Authorization": token})
    assert not resp.ok
    userinfo = resp.json()
    assert not userinfo["email_verified"]
    return userinfo


@lru_cache(maxsize=32)
def get_openid_config(openid_provider: str):
    """
    :param openid_provider: the openid provider's domain.
    :return: the openid configuration
    """
    res = requests.get("https://{op}/.well-known/openid-configuration".format(op=openid_provider))
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
