import connexion
from flask import make_response, request

from backend.common.auth0_manager import auth0_management_session
from backend.common.authorizer import get_userinfo_from_auth0
from backend.common.corpora_config import CorporaAuthConfig
from backend.common.utils.api_key import generate
from backend.common.utils.http_exceptions import NotFoundHTTPException


def get(user: str):
    print(f"\nget for key has user: {user}\n")
    print(f"\n connextion toekn info {connexion.context['token_info']}")
    print(f"\n access token maybe: {request.headers['Authorization'].split(' ')[1]}")
    identity = auth0_management_session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException()
    return make_response("", 200)


def post(user: str):
    email = get_userinfo_from_auth0(request.headers["Authorization"].split(" ")[1])["email"]
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
    api_key_id = auth0_management_session.store_api_key(password, email)
    auth0_management_session.link_api_key(user, api_key_id)
    return make_response({"key": password}, 201)


def delete(user: str):
    identity = auth0_management_session.get_user_api_key_identity(user)
    if not identity:
        raise NotFoundHTTPException()
    auth0_management_session.delete_api_key(user, identity)
    return make_response("", 202)
