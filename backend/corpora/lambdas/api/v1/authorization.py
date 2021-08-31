from functools import wraps

from flask import _request_ctx_stack
from jose import jwt

from ....common.corpora_config import CorporaAuthConfig
from .authentication import get_token
from backend.corpora.common import authorizer


# Error handler
class AuthError(Exception):
    def __init__(self, error, status_code):
        self.error = error
        self.status_code = status_code


def get_token_from_session():
    config = CorporaAuthConfig()
    try:
        token = get_token(config.cookie_name)
        return token["access_token"]
    except Exception:
        raise AuthError({"code": "token_missing", "description": "Cannot extract access token"}, 401)


def requires_auth(f):
    """Determines if the Access Token is valid"""

    @wraps(f)
    def decorated(*args, **kwargs):
        token = get_token_from_session()
        payload = authorizer.assert_authorized_token(token)

        _request_ctx_stack.top.current_user = payload
        return f(*args, **kwargs)

    return decorated


def has_scope(required_scope):
    try:
        token = get_token_from_session()
        authorizer.assert_authorized_token(token)
        unverified_claims = jwt.get_unverified_claims(token)
        if unverified_claims.get("scope"):
            token_scopes = unverified_claims["scope"].split()
            for token_scope in token_scopes:
                if token_scope == required_scope:
                    return True
            return False
    except AuthError:
        return False 


def requires_scope(required_scope):
    def actual_decorator(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            token = get_token_from_session()
            unverified_claims = jwt.get_unverified_claims(token)
            if unverified_claims.get("scope"):
                token_scopes = unverified_claims["scope"].split()
                for token_scope in token_scopes:
                    if token_scope == required_scope:
                        payload = authorizer.assert_authorized_token(token)
                        return f(user=payload["sub"], *args, **kwargs)
            raise AuthError({"code": "invalid_header", "description": "Missing permissions"}, 401)

        return wrapper

    return actual_decorator
