import json

import requests
from furl import furl
from jose import jwt, ExpiredSignatureError
from jose.exceptions import JWTClaimsError

from backend.corpora.common.authorizer import get_public_keys
from backend.corpora.common.corpora_config import CorporaAuthConfig


class Auth0APIToken:
    _bearer_token = None

    @property
    def bearer_token(self):
        if not self._bearer_token:
            self._bearer_token = self.get_bearer_token()
        else:
            unverified_header = jwt.get_unverified_header(self._bearer_token["access_token"])
            public_keys = get_public_keys(CorporaAuthConfig.api_signin_url)
            public_key = public_keys.get(unverified_header["kid"])
            try:
                jwt.decode(self._bearer_token["access_token"], public_key, algorithms=["RS256"])
            except (ExpiredSignatureError, JWTClaimsError):
                self._bearer_token = self.get_bearer_token()
        return self._bearer_token

    @staticmethod
    def get_bearer_token() -> dict:
        response = requests.post(
            CorporaAuthConfig.api_token_url,
            headers={"content-type": "application/json"},
            json={
                "client_id": CorporaAuthConfig.client_id,
                "client_secret": CorporaAuthConfig.client_secret,
                "audience": CorporaAuthConfig.api_auth0_v2_url,
                "grant_type": "client_credentials",
            },
        )
        return json.loads(response.content)["access_token"]


def get_user_info(user_id)->dict:
    """This is a functions for accessing user info from Auth0 API.
    https://auth0.com/docs/api/management/v2#!/Users/get_users_by_id

    """
    token = Auth0APIToken().bearer_token["access_token"]
    headers = {"authorization": f"Bearer {token}"}
    query_params = ",".join(["email", "given_name", "family_name"])
    path = f"api/v2/users/{user_id}"
    url = furl(CorporaAuthConfig.api_signin_url, path=path, query_params=query_params)
    response = requests.get(url.url, headers=headers)
    response.raise_for_status()
    return json.loads(response.content)
