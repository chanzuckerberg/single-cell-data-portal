from flask import make_response

from backend.corpora.common.auth0_manager import session
from backend.corpora.common.corpora_config import CorporaAuthConfig
from backend.corpora.common.utils.api_key import generate
from backend.corpora.common.utils.exceptions import NotFoundHTTPException


def get(user: str):
    identity = session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException()
    return make_response({"id": identity["username"]}, 200)


def post(user: int, token_info: dict):
    # Check if a key already exists
    identity = session.get_user_api_key_identity(user)
    if identity:
        # Delete if it exists
        session.delete_api_key(user, identity)

    # Generate a new key
    password = generate(user, CorporaAuthConfig.api_key_secret)
    key_name = password.split(".")[-1]

    api_key_id = session.store_api_key(key_name, password, token_info["email"])
    session.link_api_key(user, api_key_id)
    return make_response({"id": key_name, "api_key": password}, 202)


def delete(user: str):
    identity = session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException
    session.delete_api_key(user, identity)
    return make_response("", 201)
