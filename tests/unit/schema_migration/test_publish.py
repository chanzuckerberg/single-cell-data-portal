# ruff: noqa
from unittest import mock

from tests.unit.schema_migration.pytest_fixtures import (
    private,
    published_collection,
    revision,
    schema_migrate_and_collections,
)
from backend.layers.common.entities import DatasetProcessingStatus, DatasetStatus, DatasetVersionId


class TestPublish:
    def test_publish(self, schema_migrate_and_collections):
        schema_migrate, _ = schema_migrate_and_collections
        dataset_status = DatasetStatus.empty()
        dataset_status.processing_status = DatasetProcessingStatus.SUCCESS
        schema_migrate.business_logic.get_collection_version = mock.Mock()
        schema_migrate.business_logic.get_collection_version.return_value = mock.Mock(
            datasets=[mock.Mock(status=dataset_status)]
        )
        errors = schema_migrate.publish("collection_version_id")
        assert errors == {}
        schema_migrate.business_logic.publish_collection_version.assert_called_once()

    def test_publish__with_errors(self, schema_migrate_and_collections):
        schema_migrate, _ = schema_migrate_and_collections
        dataset_status_success = DatasetStatus.empty()
        dataset_status_success.processing_status = DatasetProcessingStatus.SUCCESS
        dataset_status_failure = DatasetStatus.empty()
        dataset_status_failure.processing_status = DatasetProcessingStatus.FAILURE
        dataset_status_failure.validation_message = "rds conversion failed"
        schema_migrate.business_logic.get_collection_version = mock.Mock()
        schema_migrate.business_logic.get_collection_version.return_value = mock.Mock(
            datasets=[
                mock.Mock(status=dataset_status_success),
                mock.Mock(version_id=DatasetVersionId("failed_dataset_version_id"), status=dataset_status_failure),
            ]
        )
        errors = schema_migrate.publish("collection_version_id")
        assert len(errors) == 1
        assert errors["failed_dataset_version_id"] == "rds conversion failed"
        schema_migrate.business_logic.publish_collection_version.assert_not_called()
