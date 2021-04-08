import unittest

import mock
import requests
from requests import HTTPError
from backend.corpora.common.utils.dl_sources.url import from_url, DropBoxURL, register, MissingHeaderException


class TestDropbox(unittest.TestCase):
    def setUp(self):
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        self.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"
        register(DropBoxURL)

    def test__get_download_url__positive(self):
        positive_tests = [
            ("https://www.dropbox.com/s/abcd1234/abcd1234?dl=0", "https://www.dropbox.com/s/abcd1234/abcd1234?dl=1"),
            ("https://www.dropbox.com/s/abcd1234/file.txt?dl=0", "https://www.dropbox.com/s/abcd1234/file.txt?dl=1"),
            ("https://www.dropbox.com/s/abcd1234/abcd1234?dl=1", "https://www.dropbox.com/s/abcd1234/abcd1234?dl=1"),
            (
                "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=1",
                "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=1",
            ),
            (
                "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx",
                "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=1",
            ),
            (
                "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=0",
                "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?dl=1",
            ),
            (
                "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?query=whatever",
                "https://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?query=whatever&dl=1",
            ),
            (
                "https://www.dropbox.com/s/abcd1234/abcd1234?dl=2",
                "https://www.dropbox.com/s/abcd1234/abcd1234?dl=2&dl=1",
            ),
            ("https://www.dropbox.com/b/abcd1234/abcd1234?dl=1", "https://www.dropbox.com/b/abcd1234/abcd1234?dl=1"),
            ("https://www.dropbox.com/s/!bcd1234/abcd1234?dl=0", "https://www.dropbox.com/s/!bcd1234/abcd1234?dl=1"),
            ("https://www.dropbox.com", "https://www.dropbox.com?dl=1"),
        ]

        for test, expected in positive_tests:
            with self.subTest(test):
                actual = from_url(test)
                self.assertEqual(expected, actual.url)

    def test_get_download_url_negative(self):
        negative_tests = [
            "https://www.notdropbox.com/s/abcd1234/abcd1234?dl=0",
            "https://www.googledrive.com/s/a1b/example.docx",
            "http://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?query=whatever&dl=1",
            "https://www.googledrive.com/s/a1b/example.docx",
            "http://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?query=whatever&dl=1",
        ]
        for test in negative_tests:
            with self.subTest(test):
                self.assertIsNone(from_url(test))

    def test__get_file_info__OK(self):
        postive_test = from_url("https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=1")
        postive_test.file_info()

    def test__get_file_info__neg(self):
        with self.subTest("404"):
            negative_test = from_url("https://www.dropbox.com/s/12345678901234/test.h5ad?dl=1")
            self.assertRaises(HTTPError, negative_test.file_info)

        with self.subTest("Missing Headers"):

            def make_response():
                response = requests.Response()
                response.status_code = 200
                response.headers = {}
                return response

            with mock.patch("requests.head", return_value=make_response()):
                negative_test = from_url("https://www.dropbox.com/s/12345678901234/test.h5ad?dl=1")
                self.assertRaises(MissingHeaderException, negative_test.file_info)
