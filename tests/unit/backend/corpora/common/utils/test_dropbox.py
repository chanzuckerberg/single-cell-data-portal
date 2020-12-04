import unittest

from requests import HTTPError

from backend.corpora.common.utils.dropbox import get_download_url_from_shared_link, get_file_info


class TestDropbox(unittest.TestCase):
    def setUp(self):
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        self.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"

    def test__get_download_url__postitive(self):
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
        ]

        for test, expected in positive_tests:
            with self.subTest(test):
                actual = get_download_url_from_shared_link(test)
                self.assertEqual(expected, actual)

    def test_get_download_url_negative(self):
        negative_tests = [
            "https://www.notdropbox.com/s/abcd1234/abcd1234?dl=0",
            "https://www.dropbox.com/b/abcd1234/abcd1234?dl=0",
            "https://www.dropbox.com/s/!bcd1234/abcd1234?dl=0",
            "https://www.dropbox.com",
            "https://www.googledrive.com/s/a1b/example.docx",
            "http://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?query=whatever&dl=1",
            "https://www.googledrive.com/s/a1b/example.docx",
            "http://www.dropbox.com/s/a1b2c3d4ef5gh6/example.docx?query=whatever&dl=1",
        ]
        for test in negative_tests:
            with self.subTest(test):
                self.assertIsNone(get_download_url_from_shared_link(test))

    def test_get_file_info(self):
        postive_test = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=1"
        get_file_info(postive_test)

        negative_test = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=1"
        self.assertRaises(HTTPError, get_file_info, negative_test)
