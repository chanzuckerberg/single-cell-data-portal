import unittest
import requests

from tests.functional.backend.corpora.common import BaseFunctionalTestCase


class TestApiKey(BaseFunctionalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

    def test_api_key_crud(self):
        headers = {"Cookie": f"cxguser={self.cookie}", "Content-Type": "application/json"}
        def _cleanup():
            requests.delete(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.addCleanup(_cleanup)

        response = requests.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)

        response = requests.post(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(201, response.status_code)
        key_1 = response.json()["key"]

        response = requests.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(200, response.status_code)

        response = requests.post(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(201, response.status_code)
        key_2 = response.json()["key"]
        self.assertNotEqual(key_1, key_2)

        response = requests.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(200, response.status_code)

        response = requests.delete(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(202, response.status_code)

        response = requests.delete(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)

        response = requests.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)


if __name__ == "__main__":
    unittest.main()
