import http.server
<<<<<<< HEAD
<<<<<<< HEAD
=======
import json
>>>>>>> Change upload to download
=======
>>>>>>> Adding dropbox specific errors
import logging
import multiprocessing
import os
import random
import socketserver
import unittest

import requests

from backend.corpora.common.corpora_orm import UploadStatus
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.math_utils import MB
from backend.corpora.dataset_processing import download


<<<<<<< HEAD
<<<<<<< HEAD
=======
logging.basicConfig(level=logging.INFO)


>>>>>>> Change upload to download
=======
>>>>>>> Adding dropbox specific errors
def start_server(path, port):
    handler = http.server.SimpleHTTPRequestHandler
    os.chdir(path)
    httpd = socketserver.TCPServer(("", port), handler)
    httpd.serve_forever()


class TestDownload(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.port = random.randint(10000, 20000)
        cls.server_process = multiprocessing.Process(
            target=start_server, args=("tests/unit/backend/corpora/fixtures", cls.port), daemon=True
        )
        cls.server_process.start()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.server_process.terminate()

    def cleanup_local_file(self, local_file):
        try:
            os.remove(local_file)
        except FileNotFoundError:
            pass

    def test_download_good(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
<<<<<<< HEAD
<<<<<<< HEAD
        url = f"http://localhost:{self.port}/upload_test_file.txt"
        file_size = int(requests.head(url).headers["content-length"])
        status = download.download("test_dataset_id", url, local_file, file_size, chunk_size=1024, update_frequency=1)
        print(status)
        self.assertEqual(1, Dataset.get("test_dataset_id").processing_status.upload_progress)
        self.assertEqual(1, status["upload_progress"])
        self.assertTrue(os.path.exists(local_file))
=======
        url = f"http://localhost:{self.port}/upload_test_file.h5ad"
=======
        url = f"http://localhost:{self.port}/upload_test_file.txt"
>>>>>>> Adding dropbox specific errors
        file_size = int(requests.head(url).headers["content-length"])
        status = download.download("test_dataset_id", url, local_file, file_size, chunk_size=1024, update_frequency=1)
        print(status)
        self.assertEqual(1, Dataset.get("test_dataset_id").processing_status.upload_progress)
        self.assertEqual(1, status["upload_progress"])
<<<<<<< HEAD
>>>>>>> Change upload to download
=======
        self.assertTrue(os.path.exists(local_file))
>>>>>>> Adding dropbox specific errors

    def test__wrong_file_size__FAILED(self):
        """Upload status is set to failed when upload progress exceeds 1. This means the file size provided is smaller
        than the file downloaded.
        """
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
<<<<<<< HEAD
<<<<<<< HEAD
        url = f"http://localhost:{self.port}/upload_test_file.txt"
=======
        url = f"http://localhost:{self.port}/upload_test_file.h5ad"
>>>>>>> Change upload to download
=======
        url = f"http://localhost:{self.port}/upload_test_file.txt"
>>>>>>> Adding dropbox specific errors

        with self.subTest("Bigger"):
            download.download("test_dataset_id", url, local_file, 1, chunk_size=1024, update_frequency=1)
            processing_status = Dataset.get("test_dataset_id").processing_status
            self.assertEqual(UploadStatus.FAILED, processing_status.upload_status)

        with self.subTest("Smaller"):
            download.download("test_dataset_id", url, local_file, 10 * MB, chunk_size=1024, update_frequency=1)
            processing_status = Dataset.get("test_dataset_id").processing_status
            self.assertEqual(UploadStatus.FAILED, processing_status.upload_status)

    def test__stop_download(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
<<<<<<< HEAD
<<<<<<< HEAD
        url = f"http://localhost:{self.port}/upload_test_file.txt"
=======
        url = f"http://localhost:{self.port}/upload_test_file.h5ad"
>>>>>>> Change upload to download
=======
        url = f"http://localhost:{self.port}/upload_test_file.txt"
>>>>>>> Adding dropbox specific errors

        progress_tracker = download.ProgressTracker(1)
        progress_tracker.stop_downloader.set()
        with self.assertLogs(download.logger, logging.INFO) as logs:
            download.downloader(url=url, local_path=local_file, chunk_size=1024, tracker=progress_tracker)
            self.assertIn("Download ended early!", logs.output[0])

    def test__bad_url__FAILED(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
<<<<<<< HEAD
<<<<<<< HEAD
        url = f"http://localhost:{self.port}/fake.txt"
=======
        url = f"http://localhost:{self.port}/fake.h5ad"
>>>>>>> Change upload to download
=======
        url = f"http://localhost:{self.port}/fake.txt"
>>>>>>> Adding dropbox specific errors
        download.download("test_dataset_id", url, local_file, 100, chunk_size=1024, update_frequency=1)
        processing_status = Dataset.get("test_dataset_id").processing_status
        self.assertEqual(UploadStatus.FAILED, processing_status.upload_status)

    def test__dataset_does_not_exist__error(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
<<<<<<< HEAD
<<<<<<< HEAD
        url = f"http://localhost:{self.port}/upload_test_file.txt"
=======
        url = f"http://localhost:{self.port}/upload_test_file.h5ad"
>>>>>>> Change upload to download
=======
        url = f"http://localhost:{self.port}/upload_test_file.txt"
>>>>>>> Adding dropbox specific errors
        file_size = int(requests.head(url).headers["content-length"])
        with self.assertRaises(AttributeError):
            download.download("test_dataset_id_fake", url, local_file, file_size, chunk_size=1024, update_frequency=1)
