import logging

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util import Retry

from backend.common.corpora_config import CorporaAuthConfig

logger = logging.getLogger(__name__)


class Auth0ManagementSession:
    """
    A wrapper around the Auth0 Management API. https://auth0.com/docs/api/management/v2

     Don't skip these tests https://github.com/chanzuckerberg/single-cell-data-portal/blob/655f391d81c89b403422ea81e0d5d69cf3427aa5/tests/functional/backend/corpora/test_api_key.py#L9
     if modifying this code.
    """

    _session = None
    _domain = None
    config = CorporaAuthConfig()

    @property
    def session(self):
        if not self._session:
            domain = self.domain
            session = requests.Session()
            retry_config = Retry(
                total=3,
                backoff_factor=1,
                status_forcelist=self.config.retry_status_forcelist,
                method_whitelist=["POST", "HEAD", "GET", "PUT", "DELETE", "OPTIONS"],
            )
            session.mount("https://", HTTPAdapter(max_retries=retry_config))

            def _refresh_token(r, *args, **kwargs):
                """Automatically refresh the auth0 management token if it expires."""
                if r.status_code == 401:
                    logger.info("Refreshing auth0 management token")
                    token = self.get_auth0_management_token(domain)
                    session.headers.update({"Authorization": token})
                    r.request.headers["Authorization"] = session.headers["Authorization"]
                    return session.send(r.request)

            session.hooks["response"].append(_refresh_token)
            self._session = session
        return self._session

    @property
    def domain(self):
        if not self._domain:
            self._domain = self.config.auth0_domain
        return self._domain

    @domain.setter
    def domain(self, domain: str):
        self._domain = domain

    def __getattr__(self, item):
        return getattr(self.session, item)

    def get_auth0_management_token(self, domain: str) -> str:
        # Generate management token
        payload = dict(
            client_id=self.config.client_id,
            client_secret=self.config.client_secret,
            grant_type="client_credentials",
            audience=f"https://{domain}/api/v2/",
        )
        headers = {
            "cache-control": "no-cache",
            "content-type": "application/json",
        }

        response = requests.post(f"https://{domain}/oauth/token", json=payload, headers=headers)
        response.raise_for_status()
        token = response.json()
        return "{} {}".format(token["token_type"], token["access_token"])

    def get_user_api_key_identity(self, primary_id: str) -> dict:
        response = self.session.get(f"https://{self.domain}/api/v2/users/{primary_id}")
        response.raise_for_status()
        body = response.json()
        identity = None
        for i in body["identities"]:
            if i["connection"] == self.config.api_key_connection_name:
                identity = i
                break
        return identity

    def delete_api_key(self, primary_id: str, identity: dict) -> None:
        user_id, provider = identity["user_id"], identity["provider"]
        # Unlink API key from user
        response = self.session.delete(
            f"https://{self.domain}/api/v2/users/{primary_id}/identities/{provider}/{user_id}"
        )
        response.raise_for_status()
        # Delete API Key
        response = self.session.delete(f"https://{self.domain}/api/v2/users/{provider}|{user_id}")
        response.raise_for_status()

    def store_api_key(self, password: str, email: str) -> str:
        # Add key to Auth0
        payload = dict(
            email=email,
            email_verified=True,
            password=password,
            connection=self.config.api_key_connection_name,
        )

        response = self.session.post(f"https://{self.domain}/api/v2/users", json=payload)
        response.raise_for_status()
        body = response.json()
        return body["identities"][0]["user_id"]

    def link_api_key(self, user: str, api_key_id: str) -> None:
        # link the user to the api key.
        payload = {"provider": "auth0", "connection_id": self.config.api_key_connection, "user_id": api_key_id}
        response = self.session.post(f"https://{self.domain}/api/v2/users/{user}/identities", json=payload)
        response.raise_for_status()

    def generate_access_token(self, user_name: str, password: str) -> dict:
        response = requests.post(
            f"https://{self.domain}/oauth/token",
            headers={"content-type": "application/x-www-form-urlencoded"},
            data=dict(
                grant_type="http://auth0.com/oauth/grant-type/password-realm",
                username=user_name,
                password=password,
                scope="profile email write:collections read:collections delete:collections",
                client_id=self.config.curator_api_client_id,
                client_secret=self.config.curator_api_client_secret,
                realm=self.config.api_key_connection_name,
                audience=self.config.curation_audience,
            ),
        )
        response.raise_for_status()
        return response.json()


auth0_management_session = Auth0ManagementSession()
