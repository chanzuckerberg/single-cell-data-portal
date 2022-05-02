from typing import Union


def has_scope(required_scope: str, scope: Union[list, str]) -> bool:
    scopes = scope.split(" ") if isinstance(scope, str) else scope
    return True if required_scope in scopes else False


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
