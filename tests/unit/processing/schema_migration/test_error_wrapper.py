import contextlib
import os
from unittest.mock import Mock


def test_error_wrapper(schema_migrate_and_collections, tmpdir):
    schema_migrate, collections = schema_migrate_and_collections
    schema_migrate.business_logic.s3_provider.upload_file = Mock()
    filename = os.path.join(tmpdir, "test_file")
    expected_file = filename + ".json"

    def func(a, b=1):
        raise Exception("error")

    with contextlib.suppress(Exception):
        schema_migrate.error_wrapper(func, filename)(1, b=2)
    assert os.path.isfile(expected_file)
    schema_migrate.business_logic.s3_provider.upload_file.assert_called_once_with(
        expected_file, "artifact-bucket", f"schema_migration/test-execution-arn/report/errors/{expected_file}", {}
    )
