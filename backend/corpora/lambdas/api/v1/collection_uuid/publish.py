from collections import defaultdict

from flask import make_response, g

from .....common.corpora_orm import CollectionVisibility
from .....common.entities import Collection, DatasetAsset
from .....common.utils.exceptions import ConflictException

from .....api_server.db import dbconnect
from .....common.utils.exceptions import ForbiddenHTTPException, ServerErrorHTTPException
from backend.corpora.lambdas.api.v1.collection import _owner_or_allowed


def check_for_duplicate_datasets(collection: Collection) -> defaultdict:
    duplicate_assets = defaultdict(list)
    etags = []
    for dataset in collection.datasets:
        if not dataset.tombstone:
            for artifact in dataset.artifacts:
                _artifact = DatasetAsset(artifact)
                metadata = _artifact.get_s3_metadata()
                if not _artifact.get_s3_metadata():
                    raise ServerErrorHTTPException("Failed to check datasets for duplications. Unable to find associated artifacts.")
                etag = metadata["ETag"]
                if etag not in etags:
                    etags.append(etag)
                else:
                    duplicate_assets[etag].append((dataset.id, _artifact.id))
    return duplicate_assets


@dbconnect
def post(collection_uuid: str, user: str):
    db_session = g.db_session
    collection = Collection.get_collection(
        db_session,
        collection_uuid,
        CollectionVisibility.PRIVATE,
        owner=_owner_or_allowed(user),
    )
    if not collection:
        raise ForbiddenHTTPException()
    if all([dataset.tombstone for dataset in collection.datasets]):
        raise ConflictException(detail="The collection must have a least one dataset.")
    duplicates = check_for_duplicate_datasets(collection)
    if duplicates:
        raise ConflictException(detail="The collection cannot have duplicate datasets.")
    collection.publish()
    return make_response({"collection_uuid": collection.id, "visibility": collection.visibility}, 202)
