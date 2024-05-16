from unittest.mock import Mock

from backend.layers.common.entities import CollectionVersionId


class TestCollectionMigrate:
    def test_can_publish_true(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        schema_migrate.schema_version = "0.0.0"
        published = collections["published"][0]
        collection_version_id = CollectionVersionId()
        schema_migrate.business_logic.create_collection_version.return_value = Mock(version_id=collection_version_id)
        datasets = [
            {
                "can_publish": "True",
                "dataset_id": dataset.dataset_id.id,
                "dataset_version_id": dataset.version_id.id,
                "collection_url": f"https://collections_domain/collections/{published.collection_id.id}",
                "execution_id": "test-execution-arn",
            }
            for dataset in published.datasets
        ]
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, True)
        assert response["collection_version_id"] == collection_version_id.id
        assert response["collection_url"] == f"https://collections_domain/collections/{published.collection_id.id}"
        assert "no_datasets" not in response
        actual_datasets = response["datasets"]
        for i in range(len(actual_datasets)):
            assert actual_datasets[i].pop("collection_version_id") == collection_version_id.id
            assert actual_datasets[i].pop("collection_id") == published.collection_id.id
            assert actual_datasets[i] == datasets[i]
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", published.collection_id.id, response
        )

    def test_can_publish_false(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        schema_migrate.schema_version = "0.0.0"
        private = collections["private"][0]
        datasets = [
            {
                "can_publish": "False",
                "dataset_id": dataset.dataset_id.id,
                "dataset_version_id": dataset.version_id.id,
                "collection_url": f"https://collections_domain/collections/{private.collection_id.id}",
                "execution_id": "test-execution-arn",
            }
            for dataset in private.datasets
        ]
        response = schema_migrate.collection_migrate(private.collection_id.id, private.version_id.id, False)
        assert response["collection_version_id"] == private.version_id.id
        assert response["collection_url"] == f"https://collections_domain/collections/{private.collection_id.id}"
        assert "no_datasets" not in response
        actual_datasets = response["datasets"]
        for i in range(len(actual_datasets)):
            assert actual_datasets[i].pop("collection_version_id") == private.version_id.id
            assert actual_datasets[i].pop("collection_id") == private.collection_id.id
            assert actual_datasets[i] == datasets[i]
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", private.collection_id.id, response
        )

    def test_can_publish_false_and_no_datasets(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        schema_migrate.schema_version = "0.0.0"
        published = collections["published"][0]
        published.datasets = []
        schema_migrate.business_logic.create_collection_version.return_value = Mock(version_id=CollectionVersionId())
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, False)
        assert response["collection_version_id"] == published.version_id.id
        assert response["collection_url"] == f"https://collections_domain/collections/{published.collection_id.id}"
        assert not response["datasets"]
        assert response["no_datasets"] == "True"
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", published.collection_id.id, response
        )

    def test_can_publish_true_and_filtered_schema_version(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        published = collections["published"][0]
        schema_migrate.business_logic.create_collection_version.return_value = Mock(version_id=CollectionVersionId())
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, False)
        assert response["collection_version_id"] == published.version_id.id
        assert response["collection_url"] == f"https://collections_domain/collections/{published.collection_id.id}"
        assert not response["datasets"]
        assert response["no_datasets"] == "True"
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", published.collection_id.id, response
        )

    def test_no_datasets(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        published = collections["published"][0]
        published.datasets = []
        schema_migrate.business_logic.create_collection_version.return_value = Mock(version_id=CollectionVersionId())
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, False)
        assert response["collection_version_id"] == published.version_id.id
        assert response["collection_url"] == f"https://collections_domain/collections/{published.collection_id.id}"
        assert not response["datasets"]
        assert response["no_datasets"] == "True"
        schema_migrate._store_sfn_response.assert_called_once_with(
            "publish_and_cleanup", published.collection_id.id, response
        )

    def test_create_migration_revision(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        schema_migrate.schema_version = "0.0.0"
        private = collections["private"][0]
        published, revision = collections["revision"]
        schema_migrate.business_logic.create_collection_version = Mock(
            return_value=Mock(version_id=CollectionVersionId())
        )

        # only call create_collection_version if the collection is published
        schema_migrate.collection_migrate(private.collection_id.id, private.version_id.id, False)
        assert not schema_migrate.business_logic.create_collection_version.called

        schema_migrate.collection_migrate(revision.collection_id.id, revision.version_id.id, False)
        assert not schema_migrate.business_logic.create_collection_version.called

        schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, False)
        schema_migrate.business_logic.create_collection_version.assert_called_once()
