import json


def mock_download_file(bucket: str, key: str, local_path: str):
    with open(local_path, "w") as f:
        json.dump({"error": key}, f)


def mock_list_directory(bucket: str, prefix: str):
    for i in range(3):
        yield f"files_{i}.json"


def test_report(schema_migrate_and_collections, tmpdir):
    schema_migrate, _ = schema_migrate_and_collections
    schema_migrate.business_logic.s3_provider.download_file = mock_download_file
    schema_migrate.business_logic.s3_provider.list_directory = mock_list_directory
    schema_migrate._upload_to_slack = lambda *args: None
    assert schema_migrate.report(tmpdir) == {
        "errors": [
            {"error": "files_0.json"},
            {"error": "files_1.json"},
            {"error": "files_2.json"},
        ]
    }
