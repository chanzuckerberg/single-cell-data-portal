import unittest

from backend.corpora.common.utils.dropbox import verify


class TestDropbox(unittest.TestCase):
    def test_regex_postitive(self):
        positive_tests = [
            "https://www.dropbox.com/s/abcd1234/abcd1234?dl=0",
            "https://www.dropbox.com/s/abcd1234/file.txt?dl=0",
            "https://www.dropbox.com/s/abcd1234/abcd1234?dl=1",
        ]

        for test in positive_tests:
            with self.subTest(test):
                self.assertTrue(verify(test))

    def test_regex_negative(self):
        negative_tests = [
            "https://www.notdropbox.com/s/abcd1234/abcd1234?dl=0",
            "https://www.dropbox.com/b/abcd1234/abcd1234?dl=0" "https://www.dropbox.com/s/!bcd1234/abcd1234?dl=0",
            "https://www.dropbox.com/s/abcd1234/abcd1234",
            "https://www.dropbox.com/s/abcd1234/abcd1234?dl=2",
            "https://www.dropbox.com",
        ]
        for test in negative_tests:
            with self.subTest(test):
                self.assertFalse(verify(test))
