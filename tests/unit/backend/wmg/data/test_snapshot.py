import os
from unittest.mock import patch

import pytest

from backend.common.census_cube.data.snapshot import (
    _get_wmg_snapshot_fullpath,
    _get_wmg_snapshot_rel_path,
    _get_wmg_snapshot_schema_dir_rel_path,
)


def test_get_wmg_snapshot_schema_dir_rel_path():
    snapshot_schema_version = "1.0.0"
    expected_path = f"snapshots/{snapshot_schema_version}"
    assert _get_wmg_snapshot_schema_dir_rel_path(snapshot_schema_version) == expected_path


def test_get_wmg_snapshot_rel_path():
    snapshot_schema_version = "1.0.0"
    snapshot_id = "test_id"
    expected_path = f"snapshots/{snapshot_schema_version}/{snapshot_id}"
    assert _get_wmg_snapshot_rel_path(snapshot_schema_version, snapshot_id) == expected_path


@pytest.mark.parametrize("snapshot_fs_root_path", [None, "/tmp"])
@patch("backend.common.census_cube.data.snapshot.CensusCubeConfig", autospec=True)
def test_get_wmg_snapshot_fullpath(mock_wmg_config, snapshot_fs_root_path):
    mock_wmg_config.return_value.bucket = "test-bucket"
    snapshot_schema_version = "1.0.0"
    snapshot_id = "test_id"
    snapshot_rel_path = f"snapshots/{snapshot_schema_version}/{snapshot_id}"
    full_path = _get_wmg_snapshot_fullpath(snapshot_rel_path, snapshot_fs_root_path)
    if snapshot_fs_root_path:
        assert full_path == os.path.join(snapshot_fs_root_path, snapshot_rel_path)
    else:
        assert full_path == os.path.join("s3://test-bucket", snapshot_rel_path)
