import pytest

from backend.layers.processing.rollback import RollbackEntity, RollbackType
from tests.unit.backend.layers.business import BaseBusinessLogicTestCase


class TestRollback:
    def test_rollback(self):
        pass

    def test_rollback__unsupported_rollback_type(self):
        pass


class TestPrivateDatasetRollback(BaseBusinessLogicTestCase):
    @pytest.mark.parametrize("pass_arg_collection_version_id", [True, False])
    def test_rollback_private_dataset(self, pass_arg_collection_version_id):
        collection_version = self.initialize_unpublished_collection(num_datasets=1)
        original_dataset_version = collection_version.datasets[0]
        new_dataset_version_id = self.business_logic.create_empty_dataset_version_for_current_dataset(
            collection_version.version_id, original_dataset_version.version_id
        ).version_id

        # Test with and without optional collection_version_id arg
        collection_version_id = collection_version.version_id if pass_arg_collection_version_id else None
        rolled_back_version = RollbackEntity(
            self.business_logic, RollbackType.PRIVATE_COLLECTIONS
        ).rollback_private_dataset(new_dataset_version_id, collection_version_id)

        # Assert returns expected rolled back version
        assert rolled_back_version.version_id.id == new_dataset_version_id.id

        # Assert DatasetVersion is rolled back
        restored_dataset_version = self.business_logic.get_collection_version(collection_version.version_id).datasets[0]
        assert restored_dataset_version.version_id.id == original_dataset_version.version_id.id

    def test_rollback_private_dataset_list(self):
        pass


class TestPrivateCollectionRollback(BaseBusinessLogicTestCase):
    def test_rollback_private_collection(self):
        pass

    def test_rollback_private_collection_list(self):
        pass


class TestPublishedCollectionRollback(BaseBusinessLogicTestCase):
    def test_rollback_public_collection(self):
        pass

    def test_rollback_published_collections(self):
        pass

    def test_rollback_published_collection_list(self):
        pass


class TestRollbackCleanup(BaseBusinessLogicTestCase):
    def test_clean_up_rolled_back_datasets(self):
        pass

    def test_clean_up_published_collection_versions(self):
        pass
