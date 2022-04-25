from typing import Union


def has_scope(required_scope: str, scopes: Union[list, str]) -> bool:
    _scopes = scopes.split(" ") if isinstance(scopes, str) else scopes
    return True if required_scope in _scopes else False


def is_user_owner_or_allowed(user: str, scope: Union[list, str], owner: str) -> bool:
    """
    Check if the user has ownership on a collection, or if it has superuser permissions
    """
    return user == owner or (has_scope("write:collections", scope))


def owner_or_allowed(user: str, scope: Union[list, str]):
    """
    Returns None if the user is superuser, `user` otherwise. Used for SQL Query where conditions
    """
    return None if has_scope("write:collections", scope) else user
