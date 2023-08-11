import unittest

import requests  # type: ignore
from jose import jwt

from tests.functional.backend.common import BaseFunctionalTestCase


class TestCollectionAccess(BaseFunctionalTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.supercurator_token = cls.get_auth_token(
            "supercurator@example.com",
            cls.config.test_auth0_user_account_password,
            additional_claims=["write:collections"],
        )
        cls.supercurator_cookie = cls.make_cookie(cls.supercurator_token)
        cls.nocollection_token = cls.get_auth_token(
            "nocollection@example.com",
            cls.config.test_auth0_user_account_password,
            additional_claims=["write:collections"],
        )
        cls.nocollection_cookie = cls.make_cookie(cls.nocollection_token)
        cls.curator_token = cls.get_auth_token(
            "nocollection@example.com",
            cls.config.test_auth0_user_account_password,
            additional_claims=["write:collections"],
        )
        cls.curator_cookie = cls.make_cookie(cls.curator_token)

    def test_collection_access(self):
        """Test that only a super curator has access to all of the collections"""
        # get collections for nocollection user
        headers = {"Cookie": f"cxguser={self.nocollection_cookie}", "Content-Type": "application/json"}
        res = self.session.get(f"{self.api}/dp/v1/collections", headers=headers)
        self.assertStatusCode(requests.codes.ok, res)
        # length should be 0
        collections = res.json()["collections"]
        private_collections = [c for c in collections if c["visibility"] == "PRIVATE"]
        self.assertEqual(len(private_collections), 0)

        # get collection for supercurator user
        headers = {"Cookie": f"cxguser={self.supercurator_cookie}", "Content-Type": "application/json"}
        res = self.session.get(f"{self.api}/dp/v1/collections", headers=headers)
        self.assertStatusCode(requests.codes.ok, res)
        # len should be a lot
        superuser_collections = [c for c in res.json()["collections"] if c["visibility"] == "PRIVATE"]

        # get collection for curator user
        headers = {"Cookie": f"cxguser={self.curator_cookie}", "Content-Type": "application/json"}
        res = self.session.get(f"{self.api}/dp/v1/collections", headers=headers)
        self.assertStatusCode(requests.codes.ok, res)

        # len should be less than super curator
        curator_collections = [c for c in res.json()["collections"] if c["visibility"] == "PRIVATE"]

        self.assertLess(len(curator_collections), len(superuser_collections))

    def test_claims(self):
        access_token = self.supercurator_token["access_token"]
        token = jwt.get_unverified_claims(access_token)
        claims = token["scope"]
        self.assertIn("write:collections", claims)

        access_token = self.curator_token["access_token"]
        token = jwt.get_unverified_claims(access_token)
        claims = token["scope"]
        self.assertNotIn("write:collections", claims)

        access_token = self.nocollection_token["access_token"]
        token = jwt.get_unverified_claims(access_token)
        claims = token["scope"]
        self.assertNotIn("write:collections", claims)


if __name__ == "__main__":
    unittest.main()
