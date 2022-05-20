from flask import make_response, jsonify

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Collection, Dataset
from backend.corpora.common.utils.http_exceptions import ForbiddenHTTPException
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed


def get_collection(db_session, collection_uuid, **kwargs):
    collection = Collection.get_collection(db_session, collection_uuid, **kwargs)
    if not collection:
        raise ForbiddenHTTPException()
    return collection


def delete_dataset_common(db_session, dataset, token_info):
    if not dataset:
        raise ForbiddenHTTPException()
    collection = Collection.get_collection(
        db_session,
        dataset.collection.id,
        owner=owner_or_allowed(token_info),
    )
    if not collection:
        raise ForbiddenHTTPException()
    if dataset.collection.visibility == CollectionVisibility.PUBLIC:
        return make_response(jsonify("Can not delete a public dataset"), 405)
    if dataset.tombstone is False:
        if dataset.published:
            dataset.update(tombstone=True, published=False)
        else:
            if dataset.original_id:
                original = Dataset.get(db_session, dataset.original_id)
                original.create_revision(dataset.collection.id)
            dataset.asset_deletion()
            dataset.delete()
