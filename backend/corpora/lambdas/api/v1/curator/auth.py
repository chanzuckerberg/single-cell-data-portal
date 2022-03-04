from flask import make_response

from backend.corpora.common.auth0_manager import session


def token_post(api_key: str):
    user_name, password = api_key.split(".")
    access_token = session.generate_access_token(user_name, password)
    return make_response(access_token, 201)
