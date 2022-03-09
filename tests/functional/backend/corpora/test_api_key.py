import unittest
import requests

from tests.functional.backend.corpora.common import BaseFunctionalTestCase


class TestApiKey(BaseFunctionalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

    def test_api_key_crud(self):
        headers = {"Cookie": f"cxguser={self.cookie}", "Content-Type": "application/json"}
        response = requests.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)

        response = requests.post(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(201, response.status_code)
        self.assertIn("key", response.json().keys())
        key_id_1 = response.json()["id"]

        response = requests.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(key_id_1, response.json()['id'])

        response = requests.post(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(201, response.status_code)
        key_id_2 = response.json()["id"]
        self.assertNotEqual(key_id_1, key_id_2)

        response = requests.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(key_id_2, response.json()['id'])

        response = requests.delete(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(202, response.status_code)

        response = requests.delete(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)

        response = requests.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertEqual(404, response.status_code)


if __name__ == "__main__":
    unittest.main()
