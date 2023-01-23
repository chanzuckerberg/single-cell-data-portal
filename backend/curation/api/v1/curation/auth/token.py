from flask import request, make_response
from jose import JWTError

from backend.common.auth0_manager import auth0_management_session
from backend.common.corpora_config import CorporaAuthConfig
from backend.common.utils import api_key
from backend.common.utils.http_exceptions import UnauthorizedError, NotFoundHTTPException


def post():
    user_api_key = request.headers["x-api-key"]
    print(f"djh user api key {user_api_key}")
    config = CorporaAuthConfig()
    try:
        token_info = api_key.verify(user_api_key, config.api_key_secret)
        print(f"djh token_info {token_info}")
    except JWTError:
        raise UnauthorizedError(detail="The API key is invalid")
    else:
        identity = auth0_management_session.get_user_api_key_identity(token_info["sub"])
        print(f"djh identity {identity}")
        if not identity:
            print(f"djh API Key is no longer valid")
            raise NotFoundHTTPException(detail="The API key is no longer valid.")
        token = auth0_management_session.generate_access_token(identity["profileData"]["email"], user_api_key)
        print(f"djh token {token}")
        return make_response(token, 201)
