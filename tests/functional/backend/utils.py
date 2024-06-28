import base64
import json
from typing import Callable, Optional

import requests
from filelock import FileLock

from tests.functional.backend.constants import AUDIENCE


def distributed_singleton(tmp_path_factory, worker_id: str, func: Callable) -> dict:
    """
    This function wraps a pytest fixture so it is only instantiated once and shared across all workers in a distributed
    test run.
    """
    if worker_id != "master":
        # not executing with multiple workers, just produce the data and let
        # pytest's fixture caching do its job
        return func()
    # get the temp directory shared by all workers
    root_tmp_dir = tmp_path_factory.getbasetemp().parent

    fn = root_tmp_dir.joinpath(func.__name__ + ".json")
    with FileLock(str(fn) + ".lock"):
        if fn.is_file():
            data = json.loads(fn.read_text())
        else:
            data = func()
            fn.write_text(json.dumps(data))
    return data


def get_auth_token(
    username: str,
    password: str,
    session: str,
    config: dict,
    deployment_stage: str,
    additional_claims: Optional[list] = None,
) -> dict[str, str]:
    standard_claims = "openid profile email offline"
    if additional_claims:
        additional_claims.append(standard_claims)
        claims = " ".join(additional_claims)
    else:
        claims = standard_claims
    response = session.post(
        "https://czi-cellxgene-dev.us.auth0.com/oauth/token",
        headers={"content-type": "application/x-www-form-urlencoded"},
        data=dict(
            grant_type="password",
            username=username,
            password=password,
            audience=AUDIENCE.get(deployment_stage),
            scope=claims,
            client_id=config.client_id,
            client_secret=config.client_secret,
        ),
    )
    response.raise_for_status()
    access_token = response.json()["access_token"]
    id_token = response.json()["id_token"]
    token = {"access_token": access_token, "id_token": id_token}
    return token


def make_cookie(auth_token: str) -> str:
    return base64.b64encode(json.dumps(auth_token).encode("utf-8")).decode()


def assertStatusCode(actual: int, expected_response: requests.Response):
    request_id = expected_response.headers.get("X-Request-Id")
    assert actual == expected_response.status_code, f"{request_id=}"
