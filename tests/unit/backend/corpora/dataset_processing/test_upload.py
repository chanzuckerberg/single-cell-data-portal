import http.server
import logging
import multiprocessing
import os
import random
import socketserver
import unittest

import requests

from backend.corpora.common.corpora_orm import UploadStatus
from backend.corpora.common.entities import Dataset
from backend.corpora.dataset_processing import upload

logging.basicConfig(level=logging.INFO)


def start_server(path, port):
    handler = http.server.SimpleHTTPRequestHandler
    os.chdir(path)
    httpd = socketserver.TCPServer(("", port), handler)
    httpd.serve_forever()


class TestUpload(unittest.TestCase):
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

    def test_upload_good(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.h5ad"
        file_size = int(requests.head(url).headers["content-length"])
        upload.upload("test_dataset_id", url, local_file, file_size, chunk_size=1024, update_frequency=1)
        self.assertTrue(os.path.exists(local_file))
        self.assertEqual(1, Dataset.get("test_dataset_id").processing_status.upload_progress)

    def test__wrong_file_size__FAILED(self):
        """Upload status is set to failed when upload progress exceeds 1. This means the file size provided is smaller
        than the file downloaded.
        """
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.h5ad"
        upload.upload("test_dataset_id", url, local_file, 100, chunk_size=1024, update_frequency=1)
        processing_status = Dataset.get("test_dataset_id").processing_status
        self.assertEqual(UploadStatus.FAILED, processing_status.upload_status)

    def test__bad_url__FAILED(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/fake.h5ad"
        upload.upload("test_dataset_id", url, local_file, 100, chunk_size=1024, update_frequency=1)
        processing_status = Dataset.get("test_dataset_id").processing_status
        self.assertEqual(UploadStatus.FAILED, processing_status.upload_status)

    def test__dataset_does_not_exist__error(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.h5ad"
        file_size = int(requests.head(url).headers["content-length"])
        with self.assertRaises(AttributeError):
            upload.upload("test_dataset_id_fake", url, local_file, file_size, chunk_size=1024, update_frequency=1)
