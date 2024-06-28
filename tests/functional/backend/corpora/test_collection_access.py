import unittest

import pytest
import requests
from jose import jwt

from tests.functional.backend.utils import assertStatusCode, distributed_singleton, get_auth_token, make_cookie


@pytest.fixture(scope="session")
def supercurator_token(session, config, deployment_stage, tmp_path_factory, worker_id):
    def _auth_token():
        return get_auth_token(
            "supercurator@example.com",
            config.test_auth0_user_account_password,
            config=config,
            session=session,
            deployment_stage=deployment_stage,
            additional_claims=["write:collections"],
        )

    return distributed_singleton(tmp_path_factory, worker_id, _auth_token)


@pytest.fixture(scope="session")
def nocollection_token(session, config, deployment_stage, tmp_path_factory, worker_id):
    def _auth_token():
        return get_auth_token(
            "nocollection@example.com",
            password=config.test_auth0_user_account_password,
            config=config,
            session=session,
            deployment_stage=deployment_stage,
            additional_claims=["write:collections"],
        )

    return distributed_singleton(tmp_path_factory, worker_id, _auth_token)


@pytest.fixture(scope="session")
def curator_token(session, config, deployment_stage, tmp_path_factory, worker_id):
    def _auth_token():
        return get_auth_token(
            "nocollection@example.com",
            config.test_auth0_user_account_password,
            config=config,
            session=session,
            deployment_stage=deployment_stage,
            additional_claims=["write:collections"],
        )

    return distributed_singleton(tmp_path_factory, worker_id, _auth_token)


@pytest.fixture(scope="session")
def supercurator_cookie(supercurator_token):
    return make_cookie(supercurator_token)


@pytest.fixture(scope="session")
def nocollection_cookie(nocollection_token):
    return make_cookie(nocollection_token)


@pytest.fixture(scope="session")
def curator_cookie(curator_token):
    return make_cookie(curator_token)


def test_collection_access(session, api_url, supercurator_cookie, nocollection_cookie, curator_cookie):
    """Test that only a super curator has access to all of the collections"""
    # get collections for nocollection user
    headers = {"Cookie": f"cxguser={nocollection_cookie}", "Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/collections", headers=headers)
    assertStatusCode(requests.codes.ok, res)
    # length should be 0
    collections = res.json()["collections"]
    private_collections = [c for c in collections if c["visibility"] == "PRIVATE"]
    assert len(private_collections) == 0

    # get collection for supercurator user
    headers = {"Cookie": f"cxguser={supercurator_cookie}", "Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/collections", headers=headers)
    assertStatusCode(requests.codes.ok, res)
    # len should be a lot
    superuser_collections = [c for c in res.json()["collections"] if c["visibility"] == "PRIVATE"]

    # get collection for curator user
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/collections", headers=headers)
    assertStatusCode(requests.codes.ok, res)

    # len should be less than super curator
    curator_collections = [c for c in res.json()["collections"] if c["visibility"] == "PRIVATE"]

    assert len(curator_collections) < len(superuser_collections)


def test_claims(supercurator_token, curator_token, nocollection_token):
    access_token = supercurator_token["access_token"]
    token = jwt.get_unverified_claims(access_token)
    claims = token["scope"]
    assert "write:collections" in claims

    access_token = curator_token["access_token"]
    token = jwt.get_unverified_claims(access_token)
    claims = token["scope"]
    assert "write:collections" not in claims

    access_token = nocollection_token["access_token"]
    token = jwt.get_unverified_claims(access_token)
    claims = token["scope"]
    assert "write:collections" not in claims


if __name__ == "__main__":
    unittest.main()
