from flask import request, make_response
from jose import JWTError

from backend.corpora.common.auth0_manager import auth0_management_session
from backend.corpora.common.corpora_config import CorporaAuthConfig
from backend.corpora.common.utils.api_key import verify
from backend.corpora.common.utils.exceptions import UnauthorizedError


def post():
    api_key = request.headers["x-api-key"]
    config = CorporaAuthConfig()
    try:
        token = verify(api_key, config.api_key_secret)
    except JWTError:
        raise UnauthorizedError("The API key is invalid")
    else:
        identity = auth0_management_session.get_user_api_key_identity(token["sub"])
        token = auth0_management_session.generate_access_token(identity["email"], api_key)
        return make_response(token, 201)
