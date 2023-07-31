from typing import Dict, List, Tuple
from unittest import mock

import pytest

from backend.layers.common.entities import (
    CollectionId,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetId,
    DatasetProcessingStatus,
    DatasetVersionId,
)
from backend.layers.processing.schema_migration import SchemaMigrate


def make_mock_dataset_version(
    dataset_id: str = None, version_id: str = None, status: dict = None, metadata: dict = None
):
    dataset_version = mock.Mock()
    dataset_version.dataset_id = DatasetId(dataset_id)
    dataset_version.version_id = DatasetVersionId(version_id)
    dataset_version.metadata.schema_version = "1.0.0"
    if metadata:
        dataset_version.metadata.configure_mock(**metadata)
    dataset_version.status.processing_status = DatasetProcessingStatus.SUCCESS
    if status:
        dataset_version.status.configure_mock(**status)
    return dataset_version


@pytest.fixture
def published_collection():
    collection = mock.Mock(spec=CollectionVersionWithDatasets, name="published")
    collection.is_published.return_value = True
    collection.is_unpublished_version.return_value = False
    collection.is_initial_unpublished_version.return_value = False
    collection.collection_id = CollectionId()
    collection.version_id = CollectionVersionId()
    collection.datasets = []
    for _i in range(2):
        collection.datasets.append(make_mock_dataset_version())
    return collection


@pytest.fixture
def revision(published_collection):
    collection = mock.Mock(spec=CollectionVersionWithDatasets, name="revision")
    collection.is_published.return_value = False
    collection.is_unpublished_version.return_value = True
    collection.is_initial_unpublished_version.return_value = False
    collection.collection_id = published_collection.collection_id
    collection.version_id = CollectionVersionId()
    collection.datasets = []
    for _dataset in published_collection.datasets:
        collection.datasets.append(make_mock_dataset_version())
    return [published_collection, collection]


@pytest.fixture
def private():
    collection = mock.Mock(spec=CollectionVersionWithDatasets, name="private")
    collection.is_published.return_value = False
    collection.is_unpublished_version.return_value = False
    collection.is_initial_unpublished_version.return_value = True
    collection.collection_id = CollectionId()
    collection.version_id = CollectionVersionId()
    collection.datasets = []
    for _i in range(2):
        collection.datasets.append(make_mock_dataset_version())
    return collection


@pytest.fixture
def schema_migrate_and_collections(published_collection, revision, private) -> Tuple[mock.Mock, Dict[str, List]]:
    business_logic = mock.Mock()

    def _get_collection_version(collection_version_id):
        db = {
            published_collection.version_id.id: published_collection,
            revision[0].version_id.id: revision[0],
            revision[1].version_id.id: revision[1],
            private.version_id.id: private,
        }
        return db[collection_version_id.id]

    business_logic.get_collection_version = _get_collection_version
    schema_migrate = SchemaMigrate(business_logic)
    return schema_migrate, {"published": [published_collection], "revision": revision, "private": [private]}
