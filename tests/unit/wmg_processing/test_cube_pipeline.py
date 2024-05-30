import contextlib
import logging
import os
import tempfile
import unittest
from unittest import mock
from unittest.mock import Mock, patch

from backend.common.census_cube.data.snapshot import _get_wmg_snapshot_schema_dir_rel_path
from backend.wmg.pipeline.constants import MAXIMUM_ADMISSIBLE_CENSUS_SCHEMA_MAJOR_VERSION
from backend.wmg.pipeline.expression_summary_and_cell_counts import create_expression_summary_and_cell_counts_cubes
from backend.wmg.pipeline.load_cube import _get_wmg_snapshot_s3_fullpath
from backend.wmg.pipeline.pipeline import logger, main


def mock_census_schema_version_unsupported(_census):
    return f"{MAXIMUM_ADMISSIBLE_CENSUS_SCHEMA_MAJOR_VERSION+1}.0.0", "2023-11-27"


@contextlib.contextmanager
def change_directory(path):
    curdir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(curdir)


class TestCubePipe(unittest.TestCase):
    @mock.patch("backend.wmg.pipeline.pipeline.notify_slack")
    @mock.patch(
        "backend.wmg.pipeline.pipeline.run_pipeline",
        new=Mock(side_effect=Exception("testing")),
    )
    def test_exception_handle_catches_errors(self, mock_notify_slack: Mock):
        with tempfile.TemporaryDirectory() as tempdir:
            last_wd = os.getcwd()
            os.chdir(tempdir)
            with self.assertLogs(logger, logging.ERROR) as logs:
                try:
                    main()
                finally:
                    os.chdir(last_wd)
            self.assertEqual(1, len(logs.records))
            mock_notify_slack.assert_called_once()

    @unittest.skip("not implemented")
    def test_pipeline_creates_files(self):
        """
        Test the pipeline successfully copies from mock s3, loads corpus and correctly builds cube/cell_ordering files
        and correctly uploads files to mock s3
        """
        raise NotImplementedError

    def test_versioned_s3_paths(self):
        """
        Tests that the path we use for writing the snapshot is the versioned path.

        NOTE: Ideally, we would want this to be a true test of the cube pipeline that actually runs the pipeline,
        mocks the S3 upload, and then verifies that the API reads are pulling from the correct mocked S3 files.
        This way, we're testing the high-level functionality of the cube pipeline and the API calls directly.
        Since our testing infra doesn't currently support S3 mocks well, the past of least resistance is to test
        the lower level functions that WMG uses to figure out where to write and read snapshots from.
        """

        wmg_bucket_name = os.environ.get("CENSUS_CUBE_BUCKET")

        # Verify that we're writing to the correct path in s3
        dest_path = _get_wmg_snapshot_s3_fullpath("v1", "snapshot-id", True)
        self.assertEqual(dest_path, f"s3://{wmg_bucket_name}/snapshots/v1/snapshot-id")

        # Verify that we're reading from the versioned path
        verioned_read_path = _get_wmg_snapshot_schema_dir_rel_path(snapshot_schema_version="v1")
        self.assertEqual(verioned_read_path, "snapshots/v1")

    def test__pipeline_fails_if_census_schema_version_unsupported(self):
        # test that the pipeline fails if the census schema version is unsupported
        with patch(
            "backend.wmg.pipeline.expression_summary_and_cell_counts.get_census_version_and_build_date",
            mock_census_schema_version_unsupported,
        ):
            try:
                create_expression_summary_and_cell_counts_cubes(corpus_path="test")
                self.fail("Expected ValueError to be raised due to unsupported census schema version.")
            except ValueError:
                pass
