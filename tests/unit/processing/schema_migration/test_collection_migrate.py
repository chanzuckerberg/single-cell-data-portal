from unittest.mock import Mock

from backend.layers.common.entities import CollectionVersionId


class TestCollectionMigrate:
    def test_migrate_published_collection(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        schema_migrate.schema_version = "0.0.0"
        published = collections["published"][0]
        collection_version_id = CollectionVersionId()
        schema_migrate.business_logic.create_collection_version.return_value = Mock(version_id=collection_version_id)
        datasets = [
            {
                "dataset_id": dataset.dataset_id.id,
                "dataset_version_id": dataset.version_id.id,
                "collection_url": f"https://collections_domain/collections/{published.collection_id.id}",
                "execution_id": "test-execution-arn",
            }
            for dataset in published.datasets
        ]
        response, response_for_log_errors_and_cleanup, response_for_span_datasets = schema_migrate.collection_migrate(
            published.collection_id.id,
            published.version_id.id,
        )
        schema_migrate._store_sfn_response.assert_any_call(
            "log_errors_and_cleanup", published.collection_id.id, response_for_log_errors_and_cleanup
        )
        schema_migrate._store_sfn_response.assert_any_call(
            "span_datasets", published.collection_id.id, response_for_span_datasets
        )
        assert response_for_log_errors_and_cleanup["collection_version_id"] == collection_version_id.id
        assert (
            response_for_log_errors_and_cleanup["collection_url"]
            == f"https://collections_domain/collections/{published.collection_id.id}"
        )
        assert "key_name" in response
        for i in range(len(response_for_span_datasets)):
            assert response_for_span_datasets[i].pop("collection_version_id") == collection_version_id.id
            assert response_for_span_datasets[i].pop("collection_id") == published.collection_id.id
            assert response_for_span_datasets[i] == datasets[i]

    def test_migrate_private_collection(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        schema_migrate.schema_version = "0.0.0"
        private = collections["private"][0]
        datasets = [
            {
                "dataset_id": dataset.dataset_id.id,
                "dataset_version_id": dataset.version_id.id,
                "collection_url": f"https://collections_domain/collections/{private.collection_id.id}",
                "execution_id": "test-execution-arn",
            }
            for dataset in private.datasets
        ]
        response, response_for_log_errors_and_cleanup, response_for_span_datasets = schema_migrate.collection_migrate(
            private.collection_id.id,
            private.version_id.id,
        )
        schema_migrate._store_sfn_response.assert_any_call(
            "log_errors_and_cleanup", private.collection_id.id, response_for_log_errors_and_cleanup
        )
        schema_migrate._store_sfn_response.assert_any_call(
            "span_datasets", private.collection_id.id, response_for_span_datasets
        )

        # verify response_for_log_errors_and_cleanup
        assert response_for_log_errors_and_cleanup["collection_version_id"] == private.version_id.id
        assert (
            response_for_log_errors_and_cleanup["collection_url"]
            == f"https://collections_domain/collections/{private.collection_id.id}"
        )

        # Verify response
        assert "key_name" in response
        assert response["collection_version_id"] == private.version_id.id
        assert response["execution_id"] == "test-execution-arn"

        # Verify response_for_span_datasets
        for i in range(len(response_for_span_datasets)):
            assert response_for_span_datasets[i].pop("collection_version_id") == private.version_id.id
            assert response_for_span_datasets[i].pop("collection_id") == private.collection_id.id
            assert response_for_span_datasets[i] == datasets[i]

    def test_filter_schema_version(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        published = collections["published"][0]
        schema_migrate.business_logic.create_collection_version.return_value = Mock(version_id=CollectionVersionId())
        response, response_for_log_errors_and_cleanup, response_for_span_datasets = schema_migrate.collection_migrate(
            published.collection_id.id,
            published.version_id.id,
        )
        schema_migrate._store_sfn_response.assert_called_once_with(
            "log_errors_and_cleanup", published.collection_id.id, response_for_log_errors_and_cleanup
        )

        # verify response_for_log_errors_and_cleanup

        assert response_for_log_errors_and_cleanup["collection_version_id"] == published.version_id.id
        assert (
            response_for_log_errors_and_cleanup["collection_url"]
            == f"https://collections_domain/collections/{published.collection_id.id}"
        )

        # verify response_for_span_datasets
        assert not response_for_span_datasets

        # verify response
        assert "key_name" not in response
        assert response["collection_version_id"] == published.version_id.id
        assert response["execution_id"] == "test-execution-arn"

    def test_no_datasets(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        published = collections["published"][0]
        published.datasets = []
        schema_migrate.business_logic.create_collection_version.return_value = Mock(version_id=CollectionVersionId())
        response, response_for_log_errors_and_cleanup, response_for_span_datasets = schema_migrate.collection_migrate(
            published.collection_id.id,
            published.version_id.id,
        )
        schema_migrate._store_sfn_response.assert_called_once_with(
            "log_errors_and_cleanup", published.collection_id.id, response_for_log_errors_and_cleanup
        )

        # verify response_for_log_errors_and_cleanup
        assert response_for_log_errors_and_cleanup["collection_version_id"] == published.version_id.id
        assert (
            response_for_log_errors_and_cleanup["collection_url"]
            == f"https://collections_domain/collections/{published.collection_id.id}"
        )

        # verify response_for_span_datasets
        assert not response_for_span_datasets

        # verify response
        assert "key_name" not in response
        assert response["collection_version_id"] == published.version_id.id
        assert response["execution_id"] == "test-execution-arn"

    def test_create_migration_revision__private(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        schema_migrate.schema_version = "0.0.0"
        private = collections["private"][0]
        schema_migrate.business_logic.create_collection_version = Mock(
            return_value=Mock(version_id=CollectionVersionId())
        )

        # only call create_collection_version if the collection is published
        schema_migrate.collection_migrate(private.collection_id.id, private.version_id.id)
        schema_migrate.business_logic.create_collection_version.assert_not_called()

    def test_create_migration_revision__published_with_revision(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        schema_migrate.schema_version = "0.0.0"
        published, revision = collections["revision"]
        schema_migrate.business_logic.create_collection_version = Mock(
            return_value=Mock(version_id=CollectionVersionId())
        )

        schema_migrate.collection_migrate(revision.collection_id.id, revision.version_id.id)
        schema_migrate.business_logic.create_collection_version.assert_not_called()

        schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id)
        schema_migrate.business_logic.create_collection_version.assert_called_once()

    def test_create_migration_revision__published_no_revision(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        schema_migrate._store_sfn_response = Mock(wraps=schema_migrate._store_sfn_response)
        schema_migrate.schema_version = "0.0.0"
        published = collections["published"][0]

        schema_migrate.business_logic.create_collection_version = Mock(
            return_value=Mock(version_id=CollectionVersionId())
        )

        schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id)
        schema_migrate.business_logic.create_collection_version.assert_called_once()
