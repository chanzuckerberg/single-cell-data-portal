import logging
import os
import pathlib
import unittest
import tempfile
from unittest import mock
from unittest.mock import Mock
import contextlib
from backend.wmg.data.cube_pipeline import main, logger
from unit.backend.wmg.fixtures.test_anndata_object import create_anndata_test_object


@contextlib.contextmanager
def change_directory(path):
    curdir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(curdir)


class TestCubePipe(unittest.TestCase):
    @mock.patch("backend.wmg.data.cube_pipeline.notify_slack")
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

    @mock.patch("backend.corpus_asset_pipelines.integrated_corpus.transform.GENE_EXPRESSION_COUNT_MIN_THRESHOLD", 1)
    def test_pipeline_creates_files(self):
        """
        Test the pipeline successfully copies from mock s3, loads corpus and correctly builds cube/cell_ordering files
        and correctly uploads files to mock s3
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            os.mkdir(f"{temp_dir}/datasets")
            for i in range(2):
                test_anndata_object = create_anndata_test_object(num_genes=5, num_cells=5)
                anndata_file_path = f"{temp_dir}/datasets/basic_test_dataset_{i}"
                os.mkdir(anndata_file_path)
                anndata_file_name = pathlib.Path(anndata_file_path, "local.h5ad")
                anndata_file_name.touch()
                test_anndata_object.write(anndata_file_name, compression="gzip")
            with change_directory(temp_dir):
                main()
