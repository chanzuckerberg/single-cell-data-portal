import contextlib
import http.server
import logging
import multiprocessing
import os
import random
import socketserver
from unittest.mock import patch

import requests

from backend.common.entities import Dataset
from backend.common.utils.math_utils import MB
from backend.portal.pipeline.processing import download
from backend.portal.pipeline.processing.exceptions import ProcessingFailed
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


def start_server(path, port):
    handler = http.server.SimpleHTTPRequestHandler
    os.chdir(path)
    httpd = socketserver.TCPServer(("", port), handler)
    httpd.serve_forever()


class TestDownload(DataPortalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.port = random.randint(10000, 20000)
        cls.server_process = multiprocessing.Process(
            target=start_server, args=("tests/unit/backend/fixtures", cls.port), daemon=True
        )
        cls.server_process.start()

    @classmethod
    def tearDownClass(cls) -> None:
        super().tearDownClass()
        cls.server_process.terminate()

    def cleanup_local_file(self, local_file):
        with contextlib.suppress(FileNotFoundError):
            os.remove(local_file)

    def test_download_good(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.txt"
        file_size = int(requests.head(url).headers["content-length"])
        status = download.download(
            "test_dataset_id",
            url,
            local_file,
            file_size,
            chunk_size=1024,
            update_frequency=1,
        )
        print(status)
        self.assertEqual(1, Dataset.get(self.session, "test_dataset_id").processing_status.upload_progress)
        self.assertEqual(1, status["upload_progress"])
        self.assertTrue(os.path.exists(local_file))

    def test__wrong_file_size__FAILED(self):
        """ProcessingFailed Exception is raised when upload progress exceeds 1, or is below 1 after upload concludes.
        This means the file size provided is smaller or larger than the file downloaded.
        """
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.txt"

        with self.subTest("Bigger"), self.assertRaises(ProcessingFailed):
            download.download(
                "test_dataset_id",
                url,
                local_file,
                1,
                chunk_size=1024,
                update_frequency=1,
            )

        with self.subTest("Smaller"), self.assertRaises(ProcessingFailed):
            download.download(
                "test_dataset_id",
                url,
                local_file,
                10 * MB,
                chunk_size=1024,
                update_frequency=1,
            )

    def test__stop_download(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.txt"

        progress_tracker = download.ProgressTracker(
            1,
        )
        progress_tracker.stop_downloader.set()
        with self.assertLogs(download.logger, logging.INFO) as logs:
            download.downloader(url=url, local_path=local_file, chunk_size=1024, tracker=progress_tracker)
            self.assertIn("Download ended early!", logs.output[1])

    def test__bad_url__FAILED(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/fake.txt"
        with self.assertRaises(ProcessingFailed):
            download.download(
                "test_dataset_id",
                url,
                local_file,
                100,
                chunk_size=1024,
                update_frequency=1,
            )

    def test__dataset_does_not_exist__error(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.txt"
        file_size = int(requests.head(url).headers["content-length"])
        with self.assertRaises(AttributeError):
            download.download(
                "test_dataset_id_fake",
                url,
                local_file,
                file_size,
                chunk_size=1024,
                update_frequency=1,
            )

    @patch("shutil.disk_usage")
    def test__disk_space__error(self, mock_disk_usage):
        mock_disk_usage.return_value = (0, 0, 0)
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.txt"
        file_size = int(requests.head(url).headers["content-length"])
        with self.assertRaises(ProcessingFailed):
            download.download("test_dataset_id_fake", url, local_file, file_size)
