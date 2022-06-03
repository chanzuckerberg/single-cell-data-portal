from sqlalchemy.orm import Session

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Collection, Dataset
from backend.corpora.common.utils.http_exceptions import ForbiddenHTTPException, MethodNotAllowedException
from backend.corpora.lambdas.api.v1.authorization import owner_or_allowed


def get_collection(db_session, collection_uuid, **kwargs):
    collection = Collection.get_collection(db_session, collection_uuid, **kwargs)
    if not collection:
        raise ForbiddenHTTPException()
    return collection


def delete_dataset_common(db_session: Session, dataset: Dataset, token_info: dict):
    if not dataset:
        raise ForbiddenHTTPException()
    collection = get_collection(
        db_session,
        dataset.collection.id,
        owner=owner_or_allowed(token_info),
    )
    if not collection:
        raise ForbiddenHTTPException()
    if dataset.collection.visibility == CollectionVisibility.PUBLIC:
        raise MethodNotAllowedException("Cannot delete a public Dataset")
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
