# ruff: noqa
import os

from unittest import mock

from tests.unit.schema_migration.pytest_fixtures import (
    private,
    published_collection,
    revision,
    schema_migrate_and_collections,
)
from backend.layers.common.entities import DatasetProcessingStatus, DatasetVersionId


@mock.patch.dict(os.environ, {"ARTIFACT_BUCKET": "upload_bucket"})
@mock.patch("backend.schema_migration.migrate.cellxgene_schema")
class TestPublishAndCleanup:
    def test_publish_and_cleanup(self, mock_cellxgene_schema, schema_migrate_and_collections):
        schema_migrate, _ = schema_migrate_and_collections
        dataset_status = mock.Mock(processing_status=DatasetProcessingStatus.SUCCESS)
        metadata = mock.Mock(schema_version="1.0.0")
        schema_migrate.business_logic.get_collection_version = mock.Mock()
        schema_migrate.business_logic.get_collection_version.return_value = mock.Mock(
            datasets=[mock.Mock(version_id=DatasetVersionId(), status=dataset_status, metadata=metadata)]
        )
        mock_cellxgene_schema.get_current_schema_version.return_value = "1.0.0"

        errors = schema_migrate.publish_and_cleanup("collection_version_id", True)
        assert errors == {}
        schema_migrate.business_logic.publish_collection_version.assert_called_once()

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
                mock.Mock(version_id=DatasetVersionId(), status=dataset_status_success, metadata=metadata_migrated),
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
        mock_cellxgene_schema.get_current_schema_version.return_value = "1.0.0"

        errors = schema_migrate.publish_and_cleanup("collection_version_id", True)
        assert len(errors) == 2
        assert errors["failed_dataset_version_id"] == "rds conversion failed"
        assert errors["non_migrated_dataset_version_id"] == "Did Not Migrate."
        schema_migrate.business_logic.publish_collection_version.assert_not_called()

    def test_publish_and_cleanup__can_not_publish(self, mock_cellxgene_schema, schema_migrate_and_collections):
        schema_migrate, _ = schema_migrate_and_collections
        dataset_status = mock.Mock(processing_status=DatasetProcessingStatus.SUCCESS)
        metadata = mock.Mock(schema_version="1.0.0")
        schema_migrate.business_logic.get_collection_version = mock.Mock()
        schema_migrate.business_logic.get_collection_version.return_value = mock.Mock(
            datasets=[mock.Mock(version_id=DatasetVersionId(), status=dataset_status, metadata=metadata)]
        )
        mock_cellxgene_schema.get_current_schema_version.return_value = "1.0.0"

        errors = schema_migrate.publish_and_cleanup("collection_version_id", False)
        assert errors == {}
        schema_migrate.business_logic.publish_collection_version.assert_not_called()
