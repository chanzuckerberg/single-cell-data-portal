"""
 Don't skip these tests https://github.com/chanzuckerberg/single-cell-data-portal/blob/82ddc2a019bb8e057bc5783e659231d4e2dc9867/tests/functional/backend/corpora/test_api_key.py#L10
 if modifying this code.
"""
from flask import make_response

from backend.common.auth0_manager import auth0_management_session
from backend.common.corpora_config import CorporaAuthConfig
from backend.common.utils.api_key import generate
from backend.common.utils.http_exceptions import NotFoundHTTPException
from backend.portal.api.app.v1.authentication import get_userinfo


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
    password = generate(user, config.api_key_secret, days_to_live)
    api_key_id = auth0_management_session.store_api_key(password, userinfo["email"])
    auth0_management_session.link_api_key(user, api_key_id)
    return make_response({"key": password}, 201)


def delete(user: str):
    identity = auth0_management_session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException()
    auth0_management_session.delete_api_key(user, identity)
    return make_response("", 202)
