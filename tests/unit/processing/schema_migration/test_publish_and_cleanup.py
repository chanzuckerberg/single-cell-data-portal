from unittest.mock import Mock, patch

from backend.layers.common.entities import CollectionId, DatasetProcessingStatus, DatasetVersionId
from tests.unit.processing.schema_migration.conftest import make_mock_dataset_version


@patch("backend.layers.processing.schema_migration.cxs_get_current_schema_version")
@patch("backend.layers.processing.schema_migration.json.dump")
class TestPublishAndCleanup:
    def test_publish_and_cleanup(self, mock_json, mock_cxs_get_current_schema_version, schema_migrate_and_collections):
        schema_migrate, _ = schema_migrate_and_collections
        mock_cxs_get_current_schema_version.return_value = "1.0.0"
        dataset_status = Mock(processing_status=DatasetProcessingStatus.SUCCESS)
        metadata = Mock(schema_version="1.0.0")
        schema_migrate.business_logic.get_collection_version = Mock()
        schema_migrate.business_logic.get_collection_version.return_value = Mock(
            datasets=[
                Mock(
                    version_id=DatasetVersionId("successful_dataset_version_id"),
                    status=dataset_status,
                    metadata=metadata,
                )
            ]
        )

        errors = schema_migrate.publish_and_cleanup("collection_version_id", True)
        assert errors == []
        schema_migrate.business_logic.publish_collection_version.assert_called_once()
        schema_migrate.s3_provider.delete_files.assert_called_once_with(
            "artifact-bucket", ["successful_dataset_version_id/migrated.h5ad"]
        )

    def test_publish_and_cleanup__with_errors(
        self, mock_json, mock_cxs_get_current_schema_version, schema_migrate_and_collections
    ):
        schema_migrate, _ = schema_migrate_and_collections
        mock_cxs_get_current_schema_version.return_value = "1.0.0"
        collection_id = CollectionId()
        failed_dataset = make_mock_dataset_version(
            version_id="failed_dataset_version_id",
            status=dict(processing_status=DatasetProcessingStatus.FAILURE, validation_message="rds conversion failed"),
            metadata=dict(schema_version="1.0.0"),
        )
        non_migrated_dataset = make_mock_dataset_version(
            version_id="non_migrated_dataset_version_id",
            metadata=dict(schema_version="0.9.0"),
        )
        datasets = [
            make_mock_dataset_version(version_id="successful_dataset_version_id"),
            failed_dataset,
            non_migrated_dataset,
        ]
        schema_migrate.business_logic.get_collection_version = Mock()
        schema_migrate.business_logic.get_collection_version.return_value = Mock(
            datasets=datasets, collection_id=collection_id
        )

        errors = schema_migrate.publish_and_cleanup("collection_version_id", True)
        assert len(errors) == 2
        assert {
            "message": "rds conversion failed",
            "dataset_processing_status": DatasetProcessingStatus.FAILURE.name,
            "collection_id": collection_id.id,
            "collection_version_id": "collection_version_id",
            "dataset_version_id": failed_dataset.version_id.id,
            "dataset_id": failed_dataset.dataset_id.id,
            "rollback": True,
        } in errors
        assert {
            "message": "Did Not Migrate.",
            "collection_id": collection_id.id,
            "collection_version_id": "collection_version_id",
            "dataset_version_id": non_migrated_dataset.version_id.id,
            "dataset_id": non_migrated_dataset.dataset_id.id,
            "rollback": False,
        } in errors
        schema_migrate.business_logic.publish_collection_version.assert_not_called()
        schema_migrate.s3_provider.delete_files.assert_called_once_with(
            "artifact-bucket",
            [
                "successful_dataset_version_id/migrated.h5ad",
                "failed_dataset_version_id/migrated.h5ad",
                "non_migrated_dataset_version_id/migrated.h5ad",
            ],
        )

    def test_publish_and_cleanup__can_not_publish(
        self, mock_json, mock_cxs_get_current_schema_version, schema_migrate_and_collections
    ):
        schema_migrate, _ = schema_migrate_and_collections
        mock_cxs_get_current_schema_version.return_value = "1.0.0"
        dataset_status = dict(processing_status=DatasetProcessingStatus.SUCCESS)
        metadata = dict(schema_version="1.0.0")
        schema_migrate.business_logic.get_collection_version = Mock()
        schema_migrate.business_logic.get_collection_version.return_value = Mock(
            datasets=[
                make_mock_dataset_version(
                    version_id="successful_dataset_version_id", status=dataset_status, metadata=metadata
                )
            ]
        )

        errors = schema_migrate.publish_and_cleanup("collection_version_id", False)
        assert errors == []
        schema_migrate.business_logic.publish_collection_version.assert_not_called()
        schema_migrate.s3_provider.delete_files.assert_called_once_with(
            "artifact-bucket", ["successful_dataset_version_id/migrated.h5ad"]
        )
