import unittest

from backend.corpora.common.utils.api_key import generate, verify


class TestAPIKey(unittest.TestCase):
    def test_generate(self):
        user = "user"
        secret = "random_secret"
        token = generate(user, secret, 1)
        verify(token, secret)


if __name__ == "__main__":
    unittest.main()
