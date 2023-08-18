from unittest import mock
from unittest.mock import patch

from backend.layers.common.entities import CollectionVersionId


@patch("backend.layers.processing.schema_migration.cxs_get_current_schema_version")
class TestCollectionMigrate:
    def test_can_publish_true(self, mock_cxs_get_current_schema_version, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = mock.Mock(wraps=schema_migrate._store_sfn_response)
        mock_cxs_get_current_schema_version.return_value = "0.0.0"
        published = collections["published"][0]
        collection_version_id = CollectionVersionId()
        schema_migrate.business_logic.create_collection_version.return_value = mock.Mock(
            version_id=collection_version_id
        )
        datasets = [
            {"can_publish": "True", "dataset_id": dataset.dataset_id.id, "dataset_version_id": dataset.version_id.id}
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
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", published.collection_id.id, response
        )

    def test_can_publish_false(self, mock_cxs_get_current_schema_version, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = mock.Mock(wraps=schema_migrate._store_sfn_response)
        mock_cxs_get_current_schema_version.return_value = "0.0.0"
        private = collections["private"][0]
        datasets = [
            {"can_publish": "False", "dataset_id": dataset.dataset_id.id, "dataset_version_id": dataset.version_id.id}
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
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", private.collection_id.id, response
        )

    def test_can_publish_false_and_no_datasets(
        self, mock_cxs_get_current_schema_version, schema_migrate_and_collections
    ):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = mock.Mock(wraps=schema_migrate._store_sfn_response)
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
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", published.collection_id.id, response
        )

    def test_can_publish_true_and_filtered_schema_version(
        self, mock_cxs_get_current_schema_version, schema_migrate_and_collections
    ):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = mock.Mock(wraps=schema_migrate._store_sfn_response)
        mock_cxs_get_current_schema_version.return_value = "1.0.0"
        published = collections["published"][0]
        schema_migrate.business_logic.create_collection_version.return_value = mock.Mock(
            version_id=CollectionVersionId()
        )
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, False)
        assert response["collection_version_id"] == published.version_id.id
        assert not response["datasets"]
        assert response["no_datasets"] == "True"
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", published.collection_id.id, response
        )

    def test_no_datasets(self, mock_cxs_get_current_schema_version, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = mock.Mock(wraps=schema_migrate._store_sfn_response)
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
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", published.collection_id.id, response
        )
