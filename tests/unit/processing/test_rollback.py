import pytest

from backend.layers.business.rollback import RollbackEntity, RollbackType
from tests.unit.backend.layers.business import BaseBusinessLogicTestCase


class TestRollback(BaseBusinessLogicTestCase):
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
