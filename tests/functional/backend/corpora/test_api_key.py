import os
import time
import unittest

import pytest
from tenacity import retry, stop_after_attempt, wait_incrementing

from tests.functional.backend.utils import assertStatusCode


@pytest.mark.skipIf(
    os.environ["DEPLOYMENT_STAGE"] == "rdev",
    "Skipping for the rdev environment to avoid a flakey race condition. Uncomment if developing this "
    "feature to run in rdev. Restore this comment before merging to main. See "
    "https://github.com/chanzuckerberg/single-cell-data-portal/issues/6198",
)
def test_api_key_crud(session, api_url, curator_cookie, request):
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}

    def _cleanup():
        session.delete(f"{api_url}/dp/v1/auth/key", headers=headers)

    request.addfinalizer(_cleanup)

    response = session.get(f"{api_url}/dp/v1/auth/key", headers=headers)
    assertStatusCode(404, response)

    response = session.post(f"{api_url}/dp/v1/auth/key", headers=headers)
    assertStatusCode(201, response)
    key_1 = response.json()["key"]

    response = session.post(
        f"{api_url}/curation/v1/auth/token",
        headers={"x-api-key": f"{key_1}", "Content-Type": "application/json"},
    )
    assertStatusCode(201, response)
    access_token = response.json()["access_token"]
    assert access_token

    # wait for auth0 User-Api-Key link to update
    @retry(wait=wait_incrementing(0, 10, 30), stop=stop_after_attempt(4))
    def get_key():
        response = session.get(f"{api_url}/dp/v1/auth/key", headers=headers)
        assertStatusCode(200, response)

    get_key()  # wait for auth0 User-Api-Key link to update

    response = session.post(f"{api_url}/dp/v1/auth/key", headers=headers)
    assertStatusCode(201, response)
    key_2 = response.json()["key"]
    assert key_1 != key_2

    # wait for auth0 User-Api-Key link to update
    time.sleep(30)

    response = session.get(f"{api_url}/dp/v1/auth/key", headers=headers)
    assertStatusCode(200, response)

    response = session.delete(f"{api_url}/dp/v1/auth/key", headers=headers)
    assertStatusCode(202, response)

    response = session.delete(f"{api_url}/dp/v1/auth/key", headers=headers)
    assertStatusCode(404, response)

    response = session.get(f"{api_url}/dp/v1/auth/key", headers=headers)
    assertStatusCode(404, response)


if __name__ == "__main__":
    unittest.main()
