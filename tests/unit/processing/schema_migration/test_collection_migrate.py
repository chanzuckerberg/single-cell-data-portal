from unittest import mock
from unittest.mock import patch

from backend.common.utils.math_utils import GB
from backend.layers.common.entities import CollectionVersionId
from backend.layers.processing.schema_migration import MAX_MEMORY, MIN_MEMORY, calculate_memory_requirements


@patch("backend.layers.processing.schema_migration.cxs_get_current_schema_version")
class TestCollectionMigrate:
    def test_can_publish_true(self, mock_cxs_get_current_schema_version, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        mock_cxs_get_current_schema_version.return_value = "0.0.0"
        published = collections["published"][0]
        collection_version_id = CollectionVersionId()
        schema_migrate.business_logic.create_collection_version.return_value = mock.Mock(
            version_id=collection_version_id
        )
        datasets = [
            {
                "can_publish": "True",
                "dataset_id": dataset.dataset_id.id,
                "dataset_version_id": dataset.version_id.id,
                "memory": "8192",
                "swap": "False",
            }
            for dataset in published.datasets
        ]
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, True)
        assert response["collection_version_id"] == collection_version_id.id
        assert "no_datasets" not in response
        actual_datasets = response["datasets"]
        for i in range(len(actual_datasets)):
            assert actual_datasets[i].pop("collection_version_id") == collection_version_id.id
            assert actual_datasets[i].pop("collection_id") == published.collection_id.id
            assert actual_datasets[i] == datasets[i]

    def test_can_publish_false(self, mock_cxs_get_current_schema_version, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        mock_cxs_get_current_schema_version.return_value = "0.0.0"
        private = collections["private"][0]
        datasets = [
            {
                "can_publish": "False",
                "dataset_id": dataset.dataset_id.id,
                "dataset_version_id": dataset.version_id.id,
                "memory": "8192",
                "swap": "False",
            }
            for dataset in private.datasets
        ]
        response = schema_migrate.collection_migrate(private.collection_id.id, private.version_id.id, False)
        assert response["collection_version_id"] == private.version_id.id
        assert "no_datasets" not in response
        actual_datasets = response["datasets"]
        for i in range(len(actual_datasets)):
            assert actual_datasets[i].pop("collection_version_id") == private.version_id.id
            assert actual_datasets[i].pop("collection_id") == private.collection_id.id
            assert actual_datasets[i] == datasets[i]

    def test_can_publish_false_and_no_datasets(
        self, mock_cxs_get_current_schema_version, schema_migrate_and_collections
    ):
        schema_migrate, collections = schema_migrate_and_collections
        mock_cxs_get_current_schema_version.return_value = "0.0.0"
        published = collections["published"][0]
        published.datasets = []
        schema_migrate.business_logic.create_collection_version.return_value = mock.Mock(
            version_id=CollectionVersionId()
        )
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, False)
        assert response["collection_version_id"] == published.version_id.id
        assert not response["datasets"]
        assert response["no_datasets"] == "True"

    def test_can_publish_true_and_filtered_schema_version(
        self, mock_cxs_get_current_schema_version, schema_migrate_and_collections
    ):
        schema_migrate, collections = schema_migrate_and_collections
        mock_cxs_get_current_schema_version.return_value = "1.0.0"
        published = collections["published"][0]
        schema_migrate.business_logic.create_collection_version.return_value = mock.Mock(
            version_id=CollectionVersionId()
        )
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, False)
        assert response["collection_version_id"] == published.version_id.id
        assert not response["datasets"]
        assert response["no_datasets"] == "True"

    def test_no_datasets(self, mock_cxs_get_current_schema_version, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        mock_cxs_get_current_schema_version.return_value = "1.0.0"
        published = collections["published"][0]
        published.datasets = []
        schema_migrate.business_logic.create_collection_version.return_value = mock.Mock(
            version_id=CollectionVersionId()
        )
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, False)
        assert response["collection_version_id"] == published.version_id.id
        assert not response["datasets"]
        assert response["no_datasets"] == "True"


class TestCalculateMemoryRequirements:
    def test_1(self):
        # Test case 1: Smaller than minimum memory (8 GB)
        dataset_size = GB  # 1 GB
        expected_memory = MIN_MEMORY * 1024  # 8 GB in MB
        expected_swap = False
        memory, swap = calculate_memory_requirements(dataset_size)
        assert memory == expected_memory
        assert swap == expected_swap

    def test_2(self):
        # Test case 2: Between minimum and maximum memory
        dataset_size = 30 * GB
        expected_memory = 30 * 2 * 1024
        expected_swap = False
        memory, swap = calculate_memory_requirements(dataset_size)
        assert memory == expected_memory
        assert swap == expected_swap

    def test_3(self):
        # Test case 3: Larger than maximum memory
        dataset_size = 100 * GB  # 100 GB
        expected_memory = MAX_MEMORY * 1024  # memory_max
        expected_swap = True
        memory, swap = calculate_memory_requirements(dataset_size)
        assert memory == expected_memory
        assert swap == expected_swap

    def test_4(self):
        # Test case 4: Dataset size is 0
        dataset_size = 0
        expected_memory = MIN_MEMORY * 1024  # Minimum memory
        expected_swap = False
        memory, swap = calculate_memory_requirements(dataset_size)
        assert memory == expected_memory
        assert swap == expected_swap
