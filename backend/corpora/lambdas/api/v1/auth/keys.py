from flask import make_response

from backend.corpora.common.auth0_manager import auth0_management_session
from backend.corpora.common.corpora_config import CorporaAuthConfig
from backend.corpora.common.utils.api_key import generate
from backend.corpora.common.utils.exceptions import NotFoundHTTPException


def get(user: str):
    identity = auth0_management_session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException()
    return make_response({"id": identity["username"]}, 200)


def post(user: str, token_info: dict):
    days_to_live = CorporaAuthConfig().days_to_live
    if not isinstance(days_to_live, (int, float)):
        days_to_live = int(days_to_live)

    # Check if a key already exists
    identity = auth0_management_session.get_user_api_key_identity(user)
    if identity:
        # Delete if it exists
        auth0_management_session.delete_api_key(user, identity)

    # Generate a new key
    password = generate(user, CorporaAuthConfig().api_key_secret, days_to_live)
    key_name = password.split(".")[-1]

    api_key_id = auth0_management_session.store_api_key(key_name, password, token_info["email"])
    auth0_management_session.link_api_key(user, api_key_id)
    return make_response({"id": key_name, "key": password}, 202)


def delete(user: str):
    identity = auth0_management_session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException
    auth0_management_session.delete_api_key(user, identity)
    return make_response("", 201)
