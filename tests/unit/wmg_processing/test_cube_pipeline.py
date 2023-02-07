import contextlib
import logging
import os
import tempfile
import unittest
from unittest import mock
from unittest.mock import Mock

from backend.wmg.pipeline.cube_pipeline import logger, main


@contextlib.contextmanager
def change_directory(path):
    curdir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(curdir)


class TestCubePipe(unittest.TestCase):
    @mock.patch("backend.wmg.pipeline.cube_pipeline.notify_slack")
    @mock.patch(
        "backend.wmg.pipeline.cube_pipeline.load_data_and_create_cube",
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
