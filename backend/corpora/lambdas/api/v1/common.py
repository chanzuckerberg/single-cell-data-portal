from sqlalchemy.orm import Session

from .authorization import is_user_owner_or_allowed
from ....common.corpora_config import CorporaConfig
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
    collection = Collection.get_collection(
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


def reshape_for_curation_api_and_is_allowed(collection, token_info, allow_access=False):
    owner = collection["owner"]
    if is_user_owner_or_allowed(token_info, owner):
        collection["access_type"] = "WRITE"
    elif not allow_access and collection["visibility"] == CollectionVisibility.PRIVATE:
        # User neither provided the uuid for access nor are they authorized by their access token
        return False
    elif token_info:
        # Access token was provided but user is not authorized
        collection["access_type"] = "READ"
    else:
        # No access token was provided
        collection["access_type"] = None

    del collection["owner"]  # Don't actually want to return 'owner' in response
    collection["collection_url"] = f"{CorporaConfig().collections_base_url}/collections/{collection['id']}"

    if "datasets" in collection:
        for dataset in collection["datasets"]:
            if "artifacts" in dataset:
                dataset["dataset_assets"] = dataset.pop("artifacts")
            if "processing_status" in dataset:
                dataset["processing_status"] = dataset["processing_status"]["processing_status"]

    return True
