from backend.corpora.common.entities import Collection
from backend.corpora.common.utils.exceptions import ForbiddenHTTPException


def get_collection(db_session, collection_uuid, **kwargs):
    collection = Collection.get_collection(db_session, collection_uuid, **kwargs)
    if not collection:
        raise ForbiddenHTTPException()
    return collection
