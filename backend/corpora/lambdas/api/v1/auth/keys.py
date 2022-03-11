from flask import make_response

from backend.corpora.common.auth0_manager import auth0_management_session
from backend.corpora.common.corpora_config import CorporaAuthConfig
from backend.corpora.common.utils.api_key import generate
from backend.corpora.common.utils.exceptions import NotFoundHTTPException
from backend.corpora.lambdas.api.v1.authentication import get_userinfo


def get(user: str):
    identity = auth0_management_session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException()
    return make_response("", 200)


def post(user: str):
    userinfo = get_userinfo()
    config = CorporaAuthConfig()
    days_to_live = config.days_to_live
    if not isinstance(days_to_live, (int, float)):
        days_to_live = int(days_to_live)

    # Check if a key already exists
    identity = auth0_management_session.get_user_api_key_identity(user)
    if identity:
        # Delete if it exists
        auth0_management_session.delete_api_key(user, identity)

    # Generate a new key
    password = generate(userinfo["email"], config().api_key_secret, days_to_live)
    api_key_id = auth0_management_session.store_api_key(password, userinfo["email"])
    auth0_management_session.link_api_key(user, api_key_id)
    return make_response({"key": password}, 202)


def delete(user: str):
    identity = auth0_management_session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException
    auth0_management_session.delete_api_key(user, identity)
    return make_response("", 201)
