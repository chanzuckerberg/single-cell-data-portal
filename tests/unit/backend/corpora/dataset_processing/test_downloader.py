import http.server
import logging
import multiprocessing
import os
import random
import requests
import socketserver

from backend.corpora.common.corpora_orm import UploadStatus
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.math_utils import MB
from backend.corpora.dataset_processing import download
from backend.corpora.dataset_processing.exceptions import ProcessingFailed
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


def start_server(path, port):
    handler = http.server.SimpleHTTPRequestHandler
    os.chdir(path)
    httpd = socketserver.TCPServer(("", port), handler)
    httpd.serve_forever()


class TestDownload(DataPortalTestCase):
    @classmethod
    def setUpClass(cls):
        DataPortalTestCase.setUpClass()
        cls.port = random.randint(10000, 20000)
        cls.server_process = multiprocessing.Process(
            target=start_server, args=("tests/unit/backend/corpora/fixtures", cls.port), daemon=True
        )
        cls.server_process.start()

    @classmethod
    def tearDownClass(cls) -> None:
        DataPortalTestCase.tearDownClass()
        cls.server_process.terminate()

    def cleanup_local_file(self, local_file):
        try:
            os.remove(local_file)
        except FileNotFoundError:
            pass

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
        """Upload status is set to failed when upload progress exceeds 1. This means the file size provided is smaller
        than the file downloaded.
        """
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.txt"

        with self.subTest("Bigger"):
            with self.assertRaises(ProcessingFailed):
                download.download(
                    "test_dataset_id",
                    url,
                    local_file,
                    1,
                    chunk_size=1024,
                    update_frequency=1,
                )
            processing_status = Dataset.get(self.session, "test_dataset_id").processing_status
            self.assertEqual(UploadStatus.FAILED, processing_status.upload_status)

        with self.subTest("Smaller"):
            with self.assertRaises(ProcessingFailed):
                download.download(
                    "test_dataset_id",
                    url,
                    local_file,
                    10 * MB,
                    chunk_size=1024,
                    update_frequency=1,
                )
            processing_status = Dataset.get(self.session, "test_dataset_id").processing_status
            self.assertEqual(UploadStatus.FAILED, processing_status.upload_status)

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
            self.assertIn("Download ended early!", logs.output[0])

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
        processing_status = Dataset.get(self.session, "test_dataset_id").processing_status
        self.assertEqual(UploadStatus.FAILED, processing_status.upload_status)

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
