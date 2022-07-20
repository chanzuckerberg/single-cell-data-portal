from typing import Union, Optional

from backend.common.utils.corpora_constants import CorporaConstants


def has_scope(required_scope: str, scope: Union[list, str]) -> bool:
    scopes = scope.split(" ") if isinstance(scope, str) else scope or []
    return True if required_scope in scopes else False


def is_super_curator(scope: Union[list, str]) -> bool:
    return has_scope(CorporaConstants.SUPER_CURATOR_SCOPE, scope)


def is_user_owner_or_allowed(user: str, scope: Union[list, str], owner: str) -> bool:
    """
    Check if the user has ownership on a collection, or if it has superuser permissions
    """
    return user == owner or is_super_curator(scope)


def owner_or_allowed(user: str, scope: Union[list, str]) -> Optional[str]:
    """
    Returns None if the user is superuser, `user` otherwise. Used for SQL Query where conditions
    """
    return None if is_super_curator(scope) else user
