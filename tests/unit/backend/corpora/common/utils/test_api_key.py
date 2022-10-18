import unittest

from backend.common.utils.api_key import generate, verify


class TestAPIKey(unittest.TestCase):
    def test_generate(self):
        user = "user"
        secret = "random_secret"
        token = generate(user, secret)
        decoded_token = verify(token, secret)
        self.assertEqual(user, decoded_token["sub"])


if __name__ == "__main__":
    unittest.main()
