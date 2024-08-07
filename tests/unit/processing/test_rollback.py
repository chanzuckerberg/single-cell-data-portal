from datetime import datetime
from unittest.mock import Mock

import pytest

from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence_mock import DatabaseProviderMock
from backend.layers.processing.rollback import RollbackEntity, RollbackType
from backend.layers.thirdparty.s3_provider_mock import MockS3Provider
from backend.layers.thirdparty.uri_provider import UriProviderInterface


@pytest.fixture
def business_logic():
    return BusinessLogic(
        DatabaseProviderMock(),
        None,
        None,
        None,
        MockS3Provider(),
        UriProviderInterface(),
    )


@pytest.fixture
def rollback_entity_private_collections(business_logic):
    return RollbackEntity(business_logic, RollbackType.PRIVATE_COLLECTIONS)


@pytest.fixture
def rollback_entity_private_collection_list(business_logic):
    return RollbackEntity(business_logic, RollbackType.PRIVATE_COLLECTION_LIST)


@pytest.fixture
def rollback_entity_public_collections(business_logic):
    return RollbackEntity(business_logic, RollbackType.PUBLIC_COLLECTIONS)


@pytest.fixture
def rollback_entity_public_collection_list(business_logic):
    return RollbackEntity(business_logic, RollbackType.PUBLIC_COLLECTION_LIST)


@pytest.fixture
def rollback_entity_private_dataset_list(business_logic):
    return RollbackEntity(business_logic, RollbackType.PRIVATE_DATASET_LIST)


def initialize_unpublished_collection(rollback_entity, num_datasets):
    version = rollback_entity.business_logic.database_provider.create_canonical_collection(
        "owner",
        "curator_name",
        None,
    )
    for _ in range(num_datasets):
        dataset_version = rollback_entity.business_logic.database_provider.create_canonical_dataset(
            version.version_id,
        )
        rollback_entity.business_logic.database_provider.add_dataset_to_collection_version_mapping(
            version.version_id, dataset_version.version_id
        )
    return rollback_entity.business_logic.database_provider.get_collection_version_with_datasets(version.version_id)


def initialize_published_collection(rollback_entity, num_datasets):
    version = initialize_unpublished_collection(rollback_entity, num_datasets=num_datasets)

    rollback_entity.business_logic.database_provider.finalize_collection_version(
        version.collection_id,
        version.version_id,
        "5.1.0",
        "1.0.0",
        published_at=datetime.utcnow(),
    )
    return rollback_entity.business_logic.database_provider.get_collection_version_with_datasets(version.version_id)


def create_and_publish_collection_revision(rollback_entity, collection_id, update_first_dataset_only=False):
    new_collection_version = rollback_entity.business_logic.create_collection_version(collection_id)

    for dataset in new_collection_version.datasets:
        rollback_entity.business_logic.create_empty_dataset_version_for_current_dataset(
            new_collection_version.version_id, dataset.version_id
        )
        if update_first_dataset_only:
            break

    rollback_entity.business_logic.database_provider.finalize_collection_version(
        collection_id,
        new_collection_version.version_id,
        "5.1.0",
        "1.0.0",
        published_at=datetime.utcnow(),
    )
    return rollback_entity.business_logic.database_provider.get_collection_version_with_datasets(
        new_collection_version.version_id
    )


# Tests


@pytest.mark.parametrize(
    "rollback_args",
    [
        ("rollback_entity_private_collections", "collections_private_rollback"),
        ("rollback_entity_private_collection_list", "collection_list_private_rollback"),
        ("rollback_entity_public_collections", "collections_public_rollback"),
        ("rollback_entity_public_collection_list", "collection_list_public_rollback"),
        ("rollback_entity_private_dataset_list", "dataset_list_private_rollback"),
    ],
)
def test_rollback(request, rollback_args):
    rollback_entity_name, rollback_function_name = rollback_args
    rollback_entity = request.getfixturevalue(rollback_entity_name)
    setattr(rollback_entity, rollback_function_name, Mock())
    rollback_entity.rollback()
    assert getattr(rollback_entity, rollback_function_name).call_count == 1


def test_rollback__unsupported_rollback_type():
    with pytest.raises(ValueError):
        RollbackEntity(Mock(), "unsupported_rollback_type").rollback()


# TestPrivateDatasetRollback


@pytest.mark.parametrize("pass_arg_collection_version_id", [True, False])
def test_rollback_private_dataset(rollback_entity_private_collections, pass_arg_collection_version_id):
    business_logic = rollback_entity_private_collections.business_logic

    collection_version = initialize_unpublished_collection(rollback_entity_private_collections, num_datasets=1)
    original_dataset_version = collection_version.datasets[0]
    newer_dataset_version_id = business_logic.create_empty_dataset_version_for_current_dataset(
        collection_version.version_id, original_dataset_version.version_id
    ).version_id
    newest_dataset_version_id = business_logic.create_empty_dataset_version_for_current_dataset(
        collection_version.version_id, newer_dataset_version_id
    ).version_id

    # Test with and without optional collection_version_id arg
    collection_version_id = collection_version.version_id if pass_arg_collection_version_id else None
    rolled_back_version = rollback_entity_private_collections.dataset_private_rollback(
        newest_dataset_version_id, collection_version_id
    )

    # Assert returns expected rolled back version
    assert rolled_back_version.version_id.id == newest_dataset_version_id.id

    # Assert DatasetVersion is rolled back to most recent previous version, not the "original"
    restored_dataset_version = business_logic.get_collection_version(collection_version.version_id).datasets[0]
    assert restored_dataset_version.version_id.id == newer_dataset_version_id.id


def test_rollback_private_dataset_list(rollback_entity_private_dataset_list):
    business_logic = rollback_entity_private_dataset_list.business_logic

    collection_version = initialize_unpublished_collection(rollback_entity_private_dataset_list, num_datasets=3)
    original_dataset_versions = collection_version.datasets
    new_dataset_version_ids = [
        business_logic.create_empty_dataset_version_for_current_dataset(
            collection_version.version_id, original_dataset_version.version_id
        ).version_id.id
        for original_dataset_version in original_dataset_versions
    ]

    # Rollback two of three new dataset versions
    rollback_entity_private_dataset_list.entity_id_list = [new_dataset_version_ids[0], new_dataset_version_ids[1]]
    rollback_entity_private_dataset_list.dataset_list_private_rollback()

    post_rollback_dataset_version_ids = [
        dataset_version.version_id.id
        for dataset_version in business_logic.get_collection_version(collection_version.version_id).datasets
    ]

    # Assert collection version is pointing to the original dataset version IDs for the rolled back datasets
    assert len(post_rollback_dataset_version_ids) == 3
    assert original_dataset_versions[0].version_id.id in post_rollback_dataset_version_ids
    assert original_dataset_versions[1].version_id.id in post_rollback_dataset_version_ids
    # Assert collection version is still pointing to newest dataset version for the non-rolled back dataset
    assert new_dataset_version_ids[2] in post_rollback_dataset_version_ids


# TestPrivateCollectionRollback


def test_rollback_private_collection(rollback_entity_private_collections):
    business_logic = rollback_entity_private_collections.business_logic

    collection_version = initialize_unpublished_collection(rollback_entity_private_collections, num_datasets=2)
    original_dataset_versions = collection_version.datasets
    newer_dataset_version_ids = [
        business_logic.create_empty_dataset_version_for_current_dataset(
            collection_version.version_id, original_dataset_version.version_id
        ).version_id
        for original_dataset_version in original_dataset_versions
    ]
    # Create a third dataset version for each dataset to roll back from, to test we rollback to "newer", not original
    for newer_dataset_version_id in newer_dataset_version_ids:
        business_logic.create_empty_dataset_version_for_current_dataset(
            collection_version.version_id, newer_dataset_version_id
        )

    rollback_entity_private_collections.collection_private_rollback(collection_version.version_id)

    post_rollback_dataset_version_ids = [
        dataset_version.version_id.id
        for dataset_version in business_logic.get_collection_version(collection_version.version_id).datasets
    ]

    # Assert collection version is pointing to the previous most recent dataset version IDs for all datasets
    assert len(post_rollback_dataset_version_ids) == 2
    assert newer_dataset_version_ids[0].id in post_rollback_dataset_version_ids
    assert newer_dataset_version_ids[1].id in post_rollback_dataset_version_ids


def test_rollback_private_collections(rollback_entity_private_collections):
    business_logic = rollback_entity_private_collections.business_logic

    original_collection_versions = [
        initialize_unpublished_collection(rollback_entity_private_collections, num_datasets=1) for _ in range(2)
    ]
    # Create a newer dataset version for each dataset to roll back from
    for collection_version in original_collection_versions:
        business_logic.create_empty_dataset_version_for_current_dataset(
            collection_version.version_id, collection_version.datasets[0].version_id
        )
    # Create published collection
    original_published_collection_version = initialize_published_collection(
        rollback_entity_private_collections, num_datasets=1
    )
    new_published_collection_version = create_and_publish_collection_revision(
        rollback_entity_private_collections, original_published_collection_version.collection_id
    )

    rollback_entity_private_collections.collections_private_rollback()

    # Assert unpublished collection versions datasets are all rolled back
    for original_collection_version in original_collection_versions:
        rolled_back_collection_version = business_logic.get_collection_version(original_collection_version.version_id)
        assert (
            rolled_back_collection_version.datasets[0].version_id.id
            == original_collection_version.datasets[0].version_id.id
        )

    # Assert published collection version is not rolled back
    published_collection = business_logic.get_canonical_collection(new_published_collection_version.collection_id)
    published_collection_version = business_logic.get_collection_version(published_collection.version_id)
    assert (
        published_collection_version.datasets[0].version_id.id
        == new_published_collection_version.datasets[0].version_id.id
    )


def test_rollback_private_collection_list(rollback_entity_private_collection_list):
    business_logic = rollback_entity_private_collection_list.business_logic

    original_collection_versions = [
        initialize_unpublished_collection(rollback_entity_private_collection_list, num_datasets=1) for _ in range(3)
    ]
    new_dataset_versions = [
        business_logic.create_empty_dataset_version_for_current_dataset(
            collection_version.version_id, collection_version.datasets[0].version_id
        ).version_id.id
        for collection_version in original_collection_versions
    ]

    rollback_entity_private_collection_list.entity_id_list = [
        original_collection_versions[0].version_id.id,
        original_collection_versions[1].version_id.id,
    ]

    rollback_entity_private_collection_list.collection_list_private_rollback()

    # Assert collection versions are pointing to the original dataset version IDs for the rolled back collections
    for original_collection_version in original_collection_versions[:2]:
        rolled_back_collection_version = business_logic.get_collection_version(original_collection_version.version_id)
        assert (
            rolled_back_collection_version.datasets[0].version_id.id
            == original_collection_version.datasets[0].version_id.id
        )

    # Assert collection version is still pointing to newest dataset version for the non-rolled back collection
    assert (
        business_logic.get_collection_version(original_collection_versions[2].version_id).datasets[0].version_id.id
        == new_dataset_versions[2]
    )


# TestPublishedCollectionRollback


def test_rollback_public_collection(rollback_entity_public_collections):
    business_logic = rollback_entity_public_collections.business_logic

    original_collection_version = initialize_published_collection(rollback_entity_public_collections, num_datasets=1)
    newer_collection_version = create_and_publish_collection_revision(
        rollback_entity_public_collections, original_collection_version.collection_id
    )
    newest_collection_version = create_and_publish_collection_revision(
        rollback_entity_public_collections, newer_collection_version.collection_id
    )

    rollback_entity_public_collections.collection_public_rollback(newest_collection_version.collection_id)

    rolled_back_collection_version = business_logic.get_canonical_collection(newest_collection_version.collection_id)

    # Assert CollectionVersionId is rolled back to most recent previous version, not the "original"
    assert rolled_back_collection_version.version_id.id == newer_collection_version.version_id.id


def test_rollback_published_collections(rollback_entity_public_collections):
    business_logic = rollback_entity_public_collections.business_logic

    original_collection_versions = [
        initialize_published_collection(rollback_entity_public_collections, num_datasets=1) for _ in range(2)
    ]

    for original_collection_version in original_collection_versions:
        create_and_publish_collection_revision(
            rollback_entity_public_collections, original_collection_version.collection_id
        )

    private_collection_version = initialize_unpublished_collection(rollback_entity_public_collections, num_datasets=1)
    new_dataset_version = business_logic.create_empty_dataset_version_for_current_dataset(
        private_collection_version.version_id, private_collection_version.datasets[0].version_id
    )

    rollback_entity_public_collections.collections_public_rollback()

    # Assert published collection versions are all rolled back
    for original_collection_version in original_collection_versions:
        rolled_back_collection = business_logic.get_canonical_collection(original_collection_version.collection_id)
        assert rolled_back_collection.version_id.id == original_collection_version.version_id.id

    # Assert private collection version is not rolled back
    private_collection = business_logic.get_collection_version(private_collection_version.version_id)
    assert private_collection.datasets[0].version_id.id == new_dataset_version.version_id.id


def test_rollback_published_collection_list(rollback_entity_public_collection_list):
    # init 2-3 public collections with 1 dataset each with prior versions. Rollback 1-2. Check those are rolled back,
    # others are pointing to new collection version still
    business_logic = rollback_entity_public_collection_list.business_logic

    original_collection_versions = [
        initialize_published_collection(rollback_entity_public_collection_list, num_datasets=1) for _ in range(3)
    ]

    new_collection_versions = [
        create_and_publish_collection_revision(
            rollback_entity_public_collection_list, original_collection_version.collection_id
        )
        for original_collection_version in original_collection_versions
    ]

    rollback_entity_public_collection_list.entity_id_list = [
        original_collection_versions[0].collection_id.id,
        original_collection_versions[1].collection_id.id,
    ]

    rollback_entity_public_collection_list.collection_list_public_rollback()

    # Assert collection versions are pointing to the original dataset version IDs for the rolled back collections
    for original_collection_version in original_collection_versions[:2]:
        rolled_back_collection = business_logic.get_canonical_collection(original_collection_version.collection_id)
        assert rolled_back_collection.version_id.id == original_collection_version.version_id.id

    # Assert collection is still pointing to newest collection version for the non-rolled back collection
    assert (
        business_logic.get_canonical_collection(original_collection_versions[2].collection_id).version_id.id
        == new_collection_versions[2].version_id.id
    )


# TestRollbackCleanUp


def test__clean_up(rollback_entity_public_collections):
    business_logic = rollback_entity_public_collections.business_logic

    original_collection_version = initialize_published_collection(rollback_entity_public_collections, num_datasets=2)

    # publish 2 new collection versions for total version history of 3
    # only update 1 of 2 datasets
    collection_revisions = [
        create_and_publish_collection_revision(
            rollback_entity_public_collections,
            original_collection_version.collection_id,
            update_first_dataset_only=True,
        )
        for _ in range(2)
    ]

    rollback_entity_public_collections.collections_public_rollback()

    rolled_back_revision = collection_revisions[-1]
    assert business_logic.get_collection_version(rolled_back_revision.version_id) is None

    revised_dataset = rolled_back_revision.datasets[0]
    unrevised_dataset = rolled_back_revision.datasets[1]
    assert business_logic.get_dataset_version(revised_dataset.version_id) is None
    assert business_logic.get_dataset_version(unrevised_dataset.version_id) is not None
