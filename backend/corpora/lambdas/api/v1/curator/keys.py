import secrets
import string

import requests
from flask import make_response

from backend.corpora.common.authorizer import CorporaAuthConfig
from backend.corpora.common.utils import json
from backend.corpora.common.utils.exceptions import NotFoundHTTPException
from backend.corpora.common.utils.singleton import Singleton


class Auth0ManagementSession(metaclass=Singleton):

    def __init__(self, domain: str = None):
        self.domain = domain if domain else CorporaAuthConfig.domain
        session = requests.Session()

        def _refresh_token(r, *args, **kwargs):
            if r.status_code == 401:
                token = Auth0ManagementSession.get_auth0_management_token()
                session.headers.update({"Authorization": token})
                r.request.headers["Authorization"] = session.headers["Authorization"]
                return session.send(r.request)

        session.hooks["response"].append(_refresh_token)
        self.session = session

    def __getattr__(self, item):
        return getattr(self, item, getattr(self.session, item))

    def get_auth0_management_token(self) -> str:
        ## Generate management token
        payload = json.dumps(
            dict(
                client_id=CorporaAuthConfig.mgmt_client_id,
                client_secret=CorporaAuthConfig.mgmt_client_secret,
                grant_type="client_credentials",
                audience=f"https://{self.domain}/api/v2/",
            )
        )
        headers = {
            "cache-control": "no-cache",
            "content-type": "application/json",
        }

        response = requests.post(f"https://{self.domain}/oauth/token", data=payload, headers=headers)
        token = response.json()
        return "{} {}".format(token["token_type"], token["access_token"])

    def get_user_api_key_identity(self, primary_id: str) -> dict:
        response = self.session.get(f"https://{domain}/api/v2/users/{primary_id}")
        response.raise_for_status()
        body = response.json()
        identity = None
        for i in body["identities"]:
            if i["connection"] == CorporaAuthConfig.api_key_connection_name:
                identity = i
        return identity

    def delete_api_key(self, primary_id: str, identity: dict) -> None:
        user_id, provider = identity["user_id"], identity["provider"]
        # Unlink API key from user
        response = self.session.delete(f"https://{domain}/api/v2/users/{primary_id}/identities/{provider}/{user_id}")
        response.raise_for_status()
        # Delete API Key
        response = self.session.delete(f"https://{domain}/api/v2/users/{provider}|{user_id}")
        response.raise_for_status()

    def generate_api_key(self, user_name: str, password: str, token_info: dict) -> str:
        # Add key to Auth0
        payload = json.dumps(
            dict(
                email=token_info["email"],
                username=user_name,
                password=password,
                connection=CorporaAuthConfig.api_key_connection_name,
            )
        )
        response = self.session.post(f"https://{domain}/api/v2/users", data=payload)
        body = response.json()
        return body["identities"][0]["user_id"]

    def link_api_key(self, user: str, api_key_id: str) -> None:
        # link the user to the api key.
        payload = json.dumps(
            {"provider": "auth0", "connection_id": CorporaAuthConfig.api_key_connection, "user_id": api_key_id}
        )
        response = self.session.post(f"https://{domain}/api/v2/users/{user}/identities", data=payload)
        response.raise_for_status()

    def generate_access_token(self, user_name: str, password: str) -> dict:
        response = self.session.post(
            f"https://{self.domain}/oauth/token",
            headers={"content-type": "application/x-www-form-urlencoded"},
            data=dict(
                grant_type="http://auth0.com/oauth/grant-type/password-realm",
                username=user_name,
                password=password,
                scope="profile email",
                client_id=CorporaAuthConfig.api_client_id,
                client_secret=CorporaAuthConfig.api_client_secret,
                realm=CorporaAuthConfig.api_key_connection,
                audience=CorporaAuthConfig.audience,
            ),
        )
        response.raise_for_status()
        return response.json()


def get(user: str):
    session = Auth0ManagementSession()
    identity = session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException()
    return make_response({"key_id": identity["username"]})


def random_string(length: int) -> str:
    return "".join(secrets.choice(string.ascii_uppercase + string.digits) for _ in range(length))


def post(user: int, token_info: dict):
    session = Auth0ManagementSession()

    # Check if a key already exists
    identity = session.get_user_api_key_identity(user)
    if identity:
        # Delete if it exists
        session.delete_api_key(user, identity)

    # Generate a new key
    user_name = random_string(8)
    password = random_string(16)

    api_key_id = session.generate_api_key(user_name, password, token_info)
    session.link_api_key(user, api_key_id)
    return make_response({"id": user_name, "api_key": f"{user_name}.{password}"}, 202)


def delete(user: str):
    session = Auth0ManagementSession()
    identity = session.get_user_api_key_identity(user)
    if not identity:
        return make_response(404)
    session.delete_api_key(user, identity)
    return make_response(201)


def token_post(api_key: str):
    session = Auth0ManagementSession()
    user_name, password = api_key.split(".")
    access_token = session.generate_access_token(user_name, password)
    return make_response(access_token, 201)
