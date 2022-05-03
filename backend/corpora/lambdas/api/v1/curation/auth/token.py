from flask import request, make_response
from jose import JWTError

from backend.corpora.common.auth0_manager import auth0_management_session
from backend.corpora.common.corpora_config import CorporaAuthConfig
from backend.corpora.common.utils import api_key
from backend.corpora.common.utils.exceptions import UnauthorizedError, NotFoundHTTPException


def post():
    user_api_key = request.headers["x-api-key"]
    config = CorporaAuthConfig()
    try:
        token_info = api_key.verify(user_api_key, config.api_key_secret)
    except JWTError:
        raise UnauthorizedError("The API key is invalid")
    else:
        identity = auth0_management_session.get_user_api_key_identity(token_info["sub"])
        if not identity:
            raise NotFoundHTTPException("The API key is no longer valid.")
        token = auth0_management_session.generate_access_token(identity["profileData"]["email"], user_api_key)
        return make_response(token, 201)
