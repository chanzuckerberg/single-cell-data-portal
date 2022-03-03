from flask import make_response

from backend.corpora.common.auth0_management_session import Auth0ManagementSession


def token_post(api_key: str):
    session = Auth0ManagementSession()
    user_name, password = api_key.split(".")
    access_token = session.generate_access_token(user_name, password)
    return make_response(access_token, 201)
