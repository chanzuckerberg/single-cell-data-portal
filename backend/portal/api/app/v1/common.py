from sqlalchemy.orm import Session

from backend.common.corpora_orm import CollectionVisibility
from backend.common.entities import Collection, Dataset
from backend.common.utils.http_exceptions import (
    ForbiddenHTTPException,
    MethodNotAllowedException,
    InvalidParametersHTTPException,
    NotFoundHTTPException,
)
from backend.portal.api.app.v1.authorization import owner_or_allowed


def get_collection_else_forbidden(db_session, collection_id, **kwargs):
    collection = Collection.get_collection(db_session, collection_id, **kwargs)
    if not collection:
        raise ForbiddenHTTPException()
    return collection


def delete_dataset_common(db_session: Session, dataset: Dataset, token_info: dict):
    if not dataset:
        raise ForbiddenHTTPException()
    get_collection_else_forbidden(
        db_session,
        dataset.collection.id,
        owner=owner_or_allowed(token_info),
    )
    if dataset.collection.visibility == CollectionVisibility.PUBLIC:
        raise MethodNotAllowedException(detail="Cannot delete a public Dataset")
    if dataset.tombstone is False:
        if dataset.published:
            dataset.update(tombstone=True, published=False)
        else:
            if dataset.original_id:
                # The dataset is a revision of a published dataset
                original = Dataset.get(db_session, dataset.original_id)
                original.create_revision(dataset.collection.id)  # Restore the original dataset and S3 assets
            dataset.asset_deletion()  # Delete the S3 assets and database rows.
            dataset.delete()  # Delete the dataset row.


def get_dataset_else_error(db_session, dataset_id, collection_id, **kwargs) -> Dataset:
    try:
        dataset = Dataset.get(db_session, dataset_id, **kwargs)
    except ValueError:
        raise InvalidParametersHTTPException()
    if not dataset:
        get_collection_else_forbidden(db_session, collection_id)  # if dataset not found, check if the collection exists
        raise NotFoundHTTPException(detail="Dataset not found.")
    return dataset
