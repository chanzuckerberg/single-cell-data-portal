from unittest import mock

from backend.layers.common.entities import DatasetProcessingStatus, DatasetVersionId


@mock.patch("backend.layers.processing.schema_migration.get_current_schema_version", return_value="1.0.0")
class TestPublishAndCleanup:
    def test_publish_and_cleanup(self, mock_get_schema_vesion, schema_migrate_and_collections):
        schema_migrate, _ = schema_migrate_and_collections
        dataset_status = mock.Mock(processing_status=DatasetProcessingStatus.SUCCESS)
        metadata = mock.Mock(schema_version="1.0.0")
        schema_migrate.business_logic.get_collection_version = mock.Mock()
        schema_migrate.business_logic.get_collection_version.return_value = mock.Mock(
            datasets=[
                mock.Mock(
                    version_id=DatasetVersionId("successful_dataset_version_id"),
                    status=dataset_status,
                    metadata=metadata,
                )
            ]
        )

        errors = schema_migrate.publish_and_cleanup("collection_version_id", True)
        assert errors == {}
        schema_migrate.business_logic.publish_collection_version.assert_called_once()
        schema_migrate.s3_provider.delete_files.assert_called_once_with(
            "artifact-bucket", ["successful_dataset_version_id/migrated.h5ad"]
        )

    def test_publish_and_cleanup__with_errors(self, mock_cellxgene_schema, schema_migrate_and_collections):
        schema_migrate, _ = schema_migrate_and_collections
        dataset_status_success = mock.Mock(processing_status=DatasetProcessingStatus.SUCCESS)
        dataset_status_failure = mock.Mock(
            processing_status=DatasetProcessingStatus.FAILURE, validation_message="rds conversion failed"
        )
        metadata_migrated = mock.Mock(schema_version="1.0.0")
        metadata_not_migrated = mock.Mock(schema_version="0.9.0")
        schema_migrate.business_logic.get_collection_version = mock.Mock()
        schema_migrate.business_logic.get_collection_version.return_value = mock.Mock(
            datasets=[
                mock.Mock(
                    version_id=DatasetVersionId("successful_dataset_version_id"),
                    status=dataset_status_success,
                    metadata=metadata_migrated,
                ),
                mock.Mock(
                    version_id=DatasetVersionId("failed_dataset_version_id"),
                    status=dataset_status_failure,
                    metadata=metadata_migrated,
                ),
                mock.Mock(
                    version_id=DatasetVersionId("non_migrated_dataset_version_id"),
                    status=dataset_status_success,
                    metadata=metadata_not_migrated,
                ),
            ]
        )

        errors = schema_migrate.publish_and_cleanup("collection_version_id", True)
        assert len(errors) == 2
        assert errors["failed_dataset_version_id"] == "rds conversion failed"
        assert errors["non_migrated_dataset_version_id"] == "Did Not Migrate."
        schema_migrate.business_logic.publish_collection_version.assert_not_called()
        schema_migrate.s3_provider.delete_files.assert_called_once_with(
            "artifact-bucket",
            [
                "successful_dataset_version_id/migrated.h5ad",
                "failed_dataset_version_id/migrated.h5ad",
                "non_migrated_dataset_version_id/migrated.h5ad",
            ],
        )

    def test_publish_and_cleanup__can_not_publish(self, mock_cellxgene_schema, schema_migrate_and_collections):
        schema_migrate, _ = schema_migrate_and_collections
        dataset_status = mock.Mock(processing_status=DatasetProcessingStatus.SUCCESS)
        metadata = mock.Mock(schema_version="1.0.0")
        schema_migrate.business_logic.get_collection_version = mock.Mock()
        schema_migrate.business_logic.get_collection_version.return_value = mock.Mock(
            datasets=[
                mock.Mock(
                    version_id=DatasetVersionId("successful_dataset_version_id"),
                    status=dataset_status,
                    metadata=metadata,
                )
            ]
        )

        errors = schema_migrate.publish_and_cleanup("collection_version_id", False)
        assert errors == {}
        schema_migrate.business_logic.publish_collection_version.assert_not_called()
        schema_migrate.s3_provider.delete_files.assert_called_once_with(
            "artifact-bucket", ["successful_dataset_version_id/migrated.h5ad"]
        )
