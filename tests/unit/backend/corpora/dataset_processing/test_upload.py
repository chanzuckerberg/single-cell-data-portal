import logging
import unittest
import http.server
import socketserver
import multiprocessing
import random
import os

from backend.corpora.common.entities import Dataset
from backend.corpora.dataset_processing import upload
from backend.corpora.common.utils.dropbox import get_file_info

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
        cls.server_thread = multiprocessing.Process(
            target=start_server, args=("tests/unit/backend/corpora/fixtures", cls.port), daemon=True
        )
        cls.server_thread.start()

    def cleanup_local_file(self, local_file):
        try:
            os.remove(local_file)
        except FileNotFoundError:
            pass

    def test_upload_good(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}/upload_test_file.h5ad"
        file_size = get_file_info(url)["size"]
        upload.upload("test_dataset_id", url, local_file, file_size, chunk_size=1024, update_frequency=1)
        self.assertTrue(os.path.exists(local_file))
        self.assertEqual(1, Dataset.get("test_dataset_id").processing_status.upload_progress)

    def test_upload_no_file_size(self):
        local_file = "local.h5ad"
        self.addCleanup(self.cleanup_local_file, local_file)
        url = f"http://localhost:{self.port}"
        file_size = -1
        upload.upload("test_dataset_id", url, local_file, file_size, chunk_size=1024, update_frequency=1)
        self.assertTrue(os.path.exists(local_file))
        self.assertEqual(-1, Dataset.get("test_dataset_id").processing_status.upload_progress)
