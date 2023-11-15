import os
import time
import unittest

from tenacity import retry, stop_after_attempt, wait_incrementing

from tests.functional.backend.common import BaseFunctionalTestCase


@unittest.skipIf(
    os.environ["DEPLOYMENT_STAGE"] == "rdev",
    "Skipping for the rdev environment to avoid a flakey race condition. Uncomment if developing this "
    "feature to run in rdev. Restore this comment before merging to main. See "
    "https://github.com/chanzuckerberg/single-cell-data-portal/issues/6198",
)
class TestApiKey(BaseFunctionalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

    def test_api_key_crud(self):
        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}

        def _cleanup():
            self.session.delete(f"{self.api}/dp/v1/auth/key", headers=headers)

        self.addCleanup(_cleanup)

        response = self.session.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertStatusCode(404, response)

        response = self.session.post(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertStatusCode(201, response)
        key_1 = response.json()["key"]

        response = self.session.post(
            f"{self.api}/curation/v1/auth/token",
            headers={"x-api-key": f"{key_1}", "Content-Type": "application/json"},
        )
        self.assertStatusCode(201, response)
        access_token = response.json()["access_token"]
        self.assertTrue(access_token)

        # wait for auth0 User-Api-Key link to update
        @retry(wait=wait_incrementing(0, 10, 30), stop=stop_after_attempt(4))
        def get_key():
            response = self.session.get(f"{self.api}/dp/v1/auth/key", headers=headers)
            self.assertStatusCode(200, response)

        get_key()  # wait for auth0 User-Api-Key link to update

        response = self.session.post(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertStatusCode(201, response)
        key_2 = response.json()["key"]
        self.assertNotEqual(key_1, key_2)

        # wait for auth0 User-Api-Key link to update
        time.sleep(30)

        response = self.session.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertStatusCode(200, response)

        response = self.session.delete(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertStatusCode(202, response)

        response = self.session.delete(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertStatusCode(404, response)

        response = self.session.get(f"{self.api}/dp/v1/auth/key", headers=headers)
        self.assertStatusCode(404, response)


if __name__ == "__main__":
    unittest.main()
