import unittest

import pytest
import requests
from jose import jwt

from tests.functional.backend.distributed import distributed_singleton
from tests.functional.backend.utils import assertStatusCode, get_auth_token, make_cookie


@pytest.fixture(scope="session")
def supercurator_token(session, config, deployment_stage, tmp_path_factory, worker_id):
    def _supercurator_token():
        return get_auth_token(
            "supercurator@example.com",
            config.test_auth0_user_account_password,
            config=config,
            session=session,
            deployment_stage=deployment_stage,
            additional_claims=["write:collections"],
        )

    return distributed_singleton(tmp_path_factory, worker_id, _supercurator_token)


@pytest.fixture(scope="session")
def nocollection_token(session, config, deployment_stage, tmp_path_factory, worker_id):
    def _nocollection_token():
        return get_auth_token(
            "nocollection@example.com",
            password=config.test_auth0_user_account_password,
            config=config,
            session=session,
            deployment_stage=deployment_stage,
            additional_claims=["write:collections"],
        )

    return distributed_singleton(tmp_path_factory, worker_id, _nocollection_token)


@pytest.fixture(scope="session")
def supercurator_cookie(supercurator_token):
    return make_cookie(supercurator_token)


@pytest.fixture(scope="session")
def nocollection_cookie(nocollection_token):
    return make_cookie(nocollection_token)


def test_nocollection_access(session, api_url, nocollection_cookie):
    """Test that a user with no private collections sees no private collections"""
    headers = {"Cookie": f"cxguser={nocollection_cookie}", "Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/collections", headers=headers)
    assertStatusCode(requests.codes.ok, res)
    collections = res.json()["collections"]
    private_collections = [c for c in collections if c["visibility"] == "PRIVATE"]
    assert len(private_collections) == 0


def test_collection_access(session, api_url, supercurator_cookie, curator_cookie):
    """Test that only a super curator has access to all of the collections"""
    # get collection for supercurator user
    headers = {"Cookie": f"cxguser={supercurator_cookie}", "Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/collections", headers=headers)
    assertStatusCode(requests.codes.ok, res)
    superuser_collections = [c for c in res.json()["collections"] if c["visibility"].lower() == "private"]

    # get collection for curator user
    headers = {"Cookie": f"cxguser={curator_cookie}", "Content-Type": "application/json"}
    res = session.get(f"{api_url}/dp/v1/collections", headers=headers)
    assertStatusCode(requests.codes.ok, res)

    curator_collections = [c for c in res.json()["collections"] if c["visibility"].lower() == "private"]

    curator_ids = {c["id"] for c in curator_collections}
    super_ids = {c["id"] for c in superuser_collections}

    print(f"[DEBUG] Super curator collection IDs ({len(super_ids)}): {sorted(super_ids)}")
    print(f"[DEBUG] Curator collection IDs ({len(curator_ids)}): {sorted(curator_ids)}")

    # Optional: explicitly print the mismatch
    if not curator_ids.issubset(super_ids):
        diff = curator_ids - super_ids
        print(f"[DEBUG] ❌ Curator sees collections that super curator cannot: {sorted(diff)}")

    assert curator_ids.issubset(
        super_ids
    ), f"Curator sees collections the super curator cannot: {curator_ids - super_ids}"
    assert len(super_ids) > len(curator_ids), "Super curator should see more collections than a curator"


def test_super_curator_claims(supercurator_token):
    access_token = supercurator_token["access_token"]
    token = jwt.get_unverified_claims(access_token)
    claims = token["scope"]
    assert "write:collections" in claims


def test_curator_claims(functest_auth_token):
    access_token = functest_auth_token["access_token"]
    token = jwt.get_unverified_claims(access_token)
    claims = token["scope"]
    assert "write:collections" not in claims


def test_nocollection_claims(nocollection_token):
    access_token = nocollection_token["access_token"]
    token = jwt.get_unverified_claims(access_token)
    claims = token["scope"]
    assert "write:collections" not in claims


if __name__ == "__main__":
    unittest.main()
