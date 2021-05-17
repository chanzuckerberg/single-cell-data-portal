import unittest

from backend.corpora.lambdas.api.v1.version import get


class TestVersion(unittest.TestCase):
    def test_get(self):
        response = get()
        self.assertIsNotNone(response)
