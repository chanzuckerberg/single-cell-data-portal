import contextlib
import logging
import os
import tempfile
import unittest
from unittest import mock
from unittest.mock import Mock

from backend.wmg.data.snapshot import _get_wmg_snapshot_schema_dir_path
from backend.wmg.pipeline import logger, main
from backend.wmg.pipeline.load_cube import _get_wmg_snapshot_s3_fullpath


@contextlib.contextmanager
def change_directory(path):
    curdir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(curdir)


class TestCubePipe(unittest.TestCase):
    @mock.patch("backend.wmg.pipeline.notify_slack")
    @mock.patch(
        "backend.wmg.pipeline.run_pipeline",
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
        Tests that the path we use for writing the snapshot is the versioned path, and verifies that the
        path we read from for API usage respects the read_versioned_snapshot param.

        NOTE: Ideally, we would want this to be a true test of the cube pipeline that actually runs the pipeline,
        mocks the S3 upload, and then verifies that the API reads are pulling from the correct mocked S3 files.
        This way, we're testing the high-level functionality of the cube pipeline and the API calls directly.
        Since our testing infra doesn't currently support S3 mocks well, the past of least resistance is to test
        the lower level functions that WMG uses to figure out where to write and read snapshots from.
        """

        wmg_bucket_name = os.environ.get("WMG_BUCKET")

        # Verify that we're writing to the correct path in s3
        dest_path = _get_wmg_snapshot_s3_fullpath("v1", "snapshot-id", True)
        self.assertEqual(dest_path, f"s3://{wmg_bucket_name}/snapshots/v1/snapshot-id")

        # Verify that we're reading from the versioned path if we pass in read_versioned_snapshot=True
        verioned_read_path = _get_wmg_snapshot_schema_dir_path(
            snapshot_schema_version="v1",
            read_versioned_snapshot=True,
        )
        self.assertEqual(verioned_read_path, "snapshots/v1")

        # Verify that we're reading from the non-versioned path if we pass in read_versioned_snapshot=False
        root_read_path = _get_wmg_snapshot_schema_dir_path(
            snapshot_schema_version="v1",
            read_versioned_snapshot=False,
        )
        self.assertEqual(root_read_path, "")
