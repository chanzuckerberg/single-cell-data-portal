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


def is_user_owner_or_allowed(token_info, owner):
    """
    Check if the user has ownership on a collection, or if it has superuser permissions
    """
    return (token_info.get("sub") and token_info.get("sub") == owner) or (has_scope("write:collections", token_info))


def owner_or_allowed(token_info):
    """
    Returns None if the user is superuser, `user` otherwise. Used for SQL Query where conditions
    """
    return None if has_scope("write:collections", token_info) else token_info.get("sub")
