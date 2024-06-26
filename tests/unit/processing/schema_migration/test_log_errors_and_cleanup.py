import json
import os
from unittest.mock import Mock, patch

from backend.layers.common.entities import DatasetProcessingStatus
from tests.unit.processing.schema_migration.conftest import make_mock_collection_version, make_mock_dataset_version


def factory_download_file():
    def download_file(bucket, key_name, local_path):
        contents = {
            "datasets": [
                {
                    "dataset_id": "dataset_id_1",
                    "dataset_version_id": "prev_successful_dataset_version_id",
                },
                {
                    "dataset_id": "dataset_id_2",
                    "dataset_version_id": "prev_failed_dataset_version_id",
                },
                {
                    "dataset_id": "dataset_id_3",
                    "dataset_version_id": "prev_non_migrated_dataset_version_id",
                },
            ],
            # these datasets populate the processed_dataset variable in the log_errors_and_cleanup function
        }
        with open(local_path, "w") as f:
            f.write(json.dumps(contents))

    return download_file


json_dump_original = json.dump


def mock_json_dump(response, f, **kwargs):
    with open(os.devnull, "w") as devnull_fp:
        json_dump_original(response, devnull_fp, **kwargs)  # tests JSON serializability of response


@patch("backend.layers.processing.schema_migration.json.dump", side_effect=mock_json_dump)
class TestLogErrorsAndCleanup:
    def test_OK(self, mock_json, schema_migrate):
        schema_migrate.business_logic.s3_provider.download_file = factory_download_file()
        schema_migrate.business_logic.s3_provider.delete_files = Mock()
        datasets = [
            make_mock_dataset_version(
                dataset_id="dataset_id_1",
                version_id="new_successful_dataset_version_id",
                status=dict(processing_status=DatasetProcessingStatus.SUCCESS),
            )
        ]
        collection_version = make_mock_collection_version(datasets)
        schema_migrate.business_logic.get_collection_version.return_value = collection_version

        errors = schema_migrate.log_errors_and_cleanup(collection_version.version_id.id)
        assert errors == []
        schema_migrate.s3_provider.delete_files.assert_any_call(
            "artifact-bucket", ["schema_migration/test-execution-arn/log_errors_and_cleanup/collection_id.json"]
        )
        schema_migrate.s3_provider.delete_files.assert_any_call(
            "artifact-bucket",
            ["prev_successful_dataset_version_id/migrated.h5ad"],
        )

    def test_with_errors(self, mock_json, schema_migrate):
        schema_migrate.business_logic.s3_provider.download_file = factory_download_file()
        schema_migrate.business_logic.s3_provider.delete_files = Mock()
        failed_dataset = make_mock_dataset_version(
            dataset_id="dataset_id_2",
            version_id="new_failed_dataset_version_id",
            status=dict(processing_status=DatasetProcessingStatus.FAILURE, validation_message="rds conversion failed"),
            metadata=dict(schema_version="1.0.0"),
        )
        non_migrated_dataset = make_mock_dataset_version(
            dataset_id="dataset_id_3",
            version_id="new_non_migrated_dataset_version_id",
            metadata=dict(schema_version="0.9.0"),
        )
        datasets = [
            failed_dataset,
            non_migrated_dataset,
        ]
        collection_version = make_mock_collection_version(datasets)
        schema_migrate.business_logic.get_collection_version.return_value = collection_version

        errors = schema_migrate.log_errors_and_cleanup(collection_version.version_id.id)
        assert len(errors) == 2
        assert {
            "message": failed_dataset.status.validation_message,
            "dataset_status": failed_dataset.status.to_dict(),
            "collection_id": collection_version.collection_id.id,
            "collection_version_id": collection_version.version_id.id,
            "dataset_version_id": failed_dataset.version_id.id,
            "dataset_id": failed_dataset.dataset_id.id,
            "rollback": True,
        } in errors
        assert {
            "message": non_migrated_dataset.status.validation_message,
            "dataset_status": non_migrated_dataset.status.to_dict(),
            "collection_id": collection_version.collection_id.id,
            "collection_version_id": collection_version.version_id.id,
            "dataset_version_id": non_migrated_dataset.version_id.id,
            "dataset_id": non_migrated_dataset.dataset_id.id,
            "rollback": False,
        } in errors
        schema_migrate.s3_provider.delete_files.assert_any_call(
            "artifact-bucket", ["schema_migration/test-execution-arn/log_errors_and_cleanup/collection_id.json"]
        )
        schema_migrate.s3_provider.delete_files.assert_any_call(
            "artifact-bucket",
            [
                "prev_failed_dataset_version_id/migrated.h5ad",
                "prev_non_migrated_dataset_version_id/migrated.h5ad",
            ],
        )

    def test_skip_unprocessed_datasets(self, mock_json, schema_migrate):
        """
        Test that datasets that do not appear in the processed_datasets variable in log_errors_and_cleanup are skipped
        """
        schema_migrate.business_logic.s3_provider.download_file = factory_download_file()
        schema_migrate.business_logic.s3_provider.delete_files = Mock()
        collection_version = make_mock_collection_version([make_mock_dataset_version()])
        schema_migrate.business_logic.get_collection_version.return_value = collection_version
        schema_migrate.check_dataset_is_latest_schema_version = Mock(return_value=True)
        errors = schema_migrate.log_errors_and_cleanup(collection_version.version_id.id)
        assert errors == []
        schema_migrate.check_dataset_is_latest_schema_version.assert_not_called()
        schema_migrate.s3_provider.delete_files.assert_any_call(
            "artifact-bucket", ["schema_migration/test-execution-arn/log_errors_and_cleanup/collection_id.json"]
        )
        schema_migrate.s3_provider.delete_files.assert_any_call("artifact-bucket", [])
