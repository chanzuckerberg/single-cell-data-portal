import datetime
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


# Tests


@pytest.mark.parametrize(
    "rollback_args",
    [
        ("rollback_entity_private_collections", "rollback_private_collections"),
        ("rollback_entity_private_collection_list", "rollback_private_collection_list"),
        ("rollback_entity_public_collections", "rollback_public_collections"),
        ("rollback_entity_public_collection_list", "rollback_public_collection_list"),
        ("rollback_entity_private_dataset_list", "rollback_private_dataset_list"),
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
    new_dataset_version_id = business_logic.create_empty_dataset_version_for_current_dataset(
        collection_version.version_id, original_dataset_version.version_id
    ).version_id

    # Test with and without optional collection_version_id arg
    collection_version_id = collection_version.version_id if pass_arg_collection_version_id else None
    rolled_back_version = rollback_entity_private_collections.rollback_private_dataset(
        new_dataset_version_id, collection_version_id
    )

    # Assert returns expected rolled back version
    assert rolled_back_version.version_id.id == new_dataset_version_id.id

    # Assert DatasetVersion is rolled back
    restored_dataset_version = business_logic.get_collection_version(collection_version.version_id).datasets[0]
    assert restored_dataset_version.version_id.id == original_dataset_version.version_id.id


def test_rollback_private_dataset_list():
    pass


# TestPrivateCollectionRollback


def test_rollback_private_collection():
    pass


def test_rollback_private_collection_list():
    pass


# TestPublishedCollectionRollback


def test_rollback_public_collection():
    pass


def test_rollback_published_collections():
    pass


def test_rollback_published_collection_list():
    pass


# TestRollbackCleanup


def test_clean_up_rolled_back_datasets():
    pass


def test_clean_up_published_collection_versions():
    pass
