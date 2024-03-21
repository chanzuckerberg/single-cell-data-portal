from typing import Dict, List, Tuple
from unittest import mock

import pytest

from backend.layers.common.entities import (
    CollectionId,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetId,
    DatasetProcessingStatus,
    DatasetStatus,
    DatasetVersionId,
)
from backend.layers.processing.schema_migration import SchemaMigrate
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProviderInterface


def make_mock_dataset_version(
    dataset_id: str = None, version_id: str = None, status: dict = None, metadata: dict = None
):
    dataset_version = mock.Mock()
    dataset_version.dataset_id = DatasetId(dataset_id)
    dataset_version.version_id = DatasetVersionId(version_id)
    dataset_version.metadata.schema_version = "1.0.0"

    # set metadata
    if metadata:
        dataset_version.metadata.configure_mock(**metadata)

    # set status
    _status = DatasetStatus.empty().to_dict()
    _status.update({"processing_status": DatasetProcessingStatus.SUCCESS} if status is None else status)
    dataset_version.status = DatasetStatus(**_status)

    return dataset_version


def make_mock_collection_version(datasets: list):
    return mock.Mock(
        datasets=datasets,
        collection_id=CollectionId("collection_id"),
        version_id=CollectionVersionId("collection_version_id"),
    )


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
def revision():
    published_collection_with_revision = mock.Mock(spec=CollectionVersionWithDatasets, name="published_with_revision")
    published_collection_with_revision.is_published.return_value = True
    published_collection_with_revision.is_unpublished_version.return_value = False
    published_collection_with_revision.is_initial_unpublished_version.return_value = False
    published_collection_with_revision.collection_id = CollectionId()
    published_collection_with_revision.version_id = CollectionVersionId()
    published_collection_with_revision.datasets = []
    for _i in range(2):
        published_collection_with_revision.datasets.append(make_mock_dataset_version())

    collection = mock.Mock(spec=CollectionVersionWithDatasets, name="revision")
    collection.is_published.return_value = False
    collection.is_unpublished_version.return_value = True
    collection.is_initial_unpublished_version.return_value = False
    collection.collection_id = published_collection_with_revision.collection_id
    collection.version_id = CollectionVersionId()
    collection.datasets = []
    for _dataset in published_collection_with_revision.datasets:
        collection.datasets.append(make_mock_dataset_version(dataset_id=_dataset.dataset_id.id))
    return [published_collection_with_revision, collection]


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
def schema_migrate(tmpdir):
    business_logic = mock.Mock()
    schema_validator = mock.Mock(spec=SchemaValidatorProviderInterface)
    schema_validator.get_current_schema_version.return_value = "1.0.0"
    schema_validator.migrate.return_value = list()
    schema_migrate = SchemaMigrate(business_logic, schema_validator)
    schema_migrate.local_path = str(tmpdir)
    return schema_migrate


@pytest.fixture
def schema_migrate_and_collections(
    tmpdir, schema_migrate, published_collection, revision, private
) -> Tuple[SchemaMigrate, Dict[str, List]]:
    db = {
        published_collection.version_id.id: published_collection,
        revision[0].version_id.id: revision[0],
        revision[1].version_id.id: revision[1],
        private.version_id.id: private,
    }

    def _get_collection_version(collection_version_id):
        return db.get(collection_version_id.id)

    def _get_collection_url(collection_id):
        return f"https://collections_domain/collections/{collection_id}"

    schema_migrate.business_logic.get_collection_version = _get_collection_version
    schema_migrate.business_logic.get_collection_url = _get_collection_url
    return schema_migrate, {"published": [published_collection], "revision": revision, "private": [private]}
