# Error handler
class AuthError(Exception):
    def __init__(self, error, status_code):
        self.error = error
        self.status_code = status_code


def has_scope(required_scope, token_info):
    try:
        if token_info.get("scope"):

            token_scopes = (
                token_info["scope"].split(" ") if isinstance(token_info["scope"], str) else token_info["scope"]
            )
            for token_scope in token_scopes:
                if token_scope == required_scope:
                    return True
        return False
    except Exception:
        return False
