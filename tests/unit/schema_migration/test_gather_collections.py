from typing import Dict, List, Tuple
from unittest import mock

import pytest

from backend.layers.common.entities import (
    CollectionId,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetVersionId,
)
from backend.schema_migration.migrate import gather_collections


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
        dataset_version = mock.Mock()
        dataset_version.dataset_id = DatasetVersionId()
        collection.datasets.append(dataset_version)
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
    for dataset in published_collection.datasets:
        dataset_version = mock.Mock()
        dataset_version.dataset_id = dataset.dataset_id
        collection.datasets.append(dataset_version)
    return [published_collection, collection]


@pytest.fixture
def private(published_collection):
    collection = mock.Mock(spec=CollectionVersionWithDatasets, name="private")
    collection.is_published.return_value = False
    collection.is_unpublished_version.return_value = False
    collection.is_initial_unpublished_version.return_value = True
    collection.collection_id = CollectionId()
    collection.version_id = CollectionVersionId()
    collection.datasets = []
    for _i in range(2):
        dataset_version = mock.Mock()
        dataset_version.dataset_id = DatasetVersionId()
        collection.datasets.append(dataset_version)
    return collection


@pytest.fixture
def business_logic_and_collections(published_collection, revision, private) -> Tuple[mock.Mock, Dict[str, List]]:
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
    return business_logic, {"published": [published_collection], "revision": revision, "private": [private]}


def test_with_revision(business_logic_and_collections):
    business_logic, collections = business_logic_and_collections
    published, revision = collections["revision"]
    with mock.patch("backend.schema_migration.migrate.get_business_logic", return_value=business_logic):
        business_logic.get_collections.side_effect = [[published], [revision]]  # get_collectioins is called twice
        response = gather_collections()
        assert published.collection_id.id not in response["published"]
        assert revision.version_id.id in response["revision"]
        assert [ds.dataset_id.id for ds in revision.datasets] == response["revision"][revision.version_id.id]


def test_with_published(business_logic_and_collections):
    business_logic, collections = business_logic_and_collections
    published = collections["published"]
    with mock.patch("backend.schema_migration.migrate.get_business_logic", return_value=business_logic):
        business_logic.get_collections.side_effect = [published, []]  # get_collectioins is called twice
        response = gather_collections()
        assert published[0].collection_id.id in response["published"]
        assert [ds.dataset_id.id for ds in published[0].datasets] == response["published"][
            published[0].collection_id.id
        ]


def test_with_private(business_logic_and_collections):
    business_logic, collections = business_logic_and_collections
    private = collections["private"]
    with mock.patch("backend.schema_migration.migrate.get_business_logic", return_value=business_logic):
        business_logic.get_collections.side_effect = [[], private]  # get_collectioins is called twice
        response = gather_collections()
        assert private[0].collection_id.id in response["private"]
        assert [ds.dataset_id.id for ds in private[0].datasets] == response["private"][private[0].collection_id.id]
