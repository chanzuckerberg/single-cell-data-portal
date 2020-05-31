import unittest

from moto import mock_s3

from backend.corpora.common.utils.s3_utils import generate_file_url


class TestS3Utils(unittest.TestCase):
    def setUp(self):
        self.s3_mock = mock_s3()
        self.s3_mock.start()

    def tearDown(self):
        self.s3_mock.stop()

    def test_generate_file_url(self):
        file_prefix = "test_prefix"

        with self.subTest("Generate URL OK"):
            url = generate_file_url(file_prefix)
            self.assertIn(file_prefix, url)
            self.assertIn("Expires=", url)
