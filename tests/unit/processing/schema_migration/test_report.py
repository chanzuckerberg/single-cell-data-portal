import json


def mock_download_file(bucket: str, key: str, local_path: str):
    with open(local_path, "w") as f:
        json.dump(key, f)


def mock_list_directory(bucket: str, prefix: str):
    if "errors" in prefix:
        for i in range(3):
            yield f"files_{i}.json"
    elif "migrate_changes" in prefix:
        for i in range(2):
            yield f"dataset_{i}_changes.json"


def test_report(schema_migrate_and_collections, tmpdir):
    schema_migrate, _ = schema_migrate_and_collections
    schema_migrate.business_logic.s3_provider.download_file = mock_download_file
    schema_migrate.business_logic.s3_provider.list_directory = mock_list_directory
    schema_migrate._upload_to_slack = lambda *args: None
    schema_migrate.local_path = str(tmpdir)
    assert schema_migrate.report(dry_run=False) == {
        "errors": ["files_0.json", "files_1.json", "files_2.json"],
        "migrate_changes": ["dataset_0_changes.json", "dataset_1_changes.json"],
    }
    schema_migrate.s3_provider.delete_prefix.assert_called_once_with(
        schema_migrate.artifact_bucket, f"schema_migration/{schema_migrate.execution_id}"
    )
