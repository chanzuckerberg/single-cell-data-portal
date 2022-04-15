from backend.corpora.common.entities import Collection
from backend.corpora.common.utils.exceptions import ForbiddenHTTPException
from backend.corpora.lambdas.api.v1.authorization import has_scope


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


def get_collection(db_session, collection_uuid, **kwargs):
    collection = Collection.get_collection(db_session, collection_uuid, **kwargs)
    if not collection:
        raise ForbiddenHTTPException()
    return collection
