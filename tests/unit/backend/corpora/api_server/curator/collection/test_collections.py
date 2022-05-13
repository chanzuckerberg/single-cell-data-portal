import unittest
from unittest.mock import patch, Mock

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.lambdas.api.v1.authentication import decode_token
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token, make_token


class TestAuthToken(BaseAuthAPITest):
    @patch("backend.corpora.lambdas.api.v1.curation.collection.dataset.sts_client")
    def test__generate_s3_credentials__OK(self, sts_client: Mock):
        def _test(token_claims: dict, additional_scope: list = None):
            token = make_token(token_claims, additional_scope=additional_scope, token_duration=10)
            sts_client.assume_role_with_web_identity = Mock(
                return_value={
                    "access_key": "test_key",
                    "secret_access_key": "test_session_token",
                    "session_token": "test_session_token",
                }
            )
            collection = self.generate_collection(self.session)
            headers = {"Authorization": f"Bearer {token}"}

            response = self.app.post(
                f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials", headers=headers
            )
            self.assertEqual(200, response.status_code)

        with self.subTest("collection owner"):
            _test(
                dict(
                    sub="test_user_id",
                    email="fake_user@email.com",
                )
            )

        with self.subTest("super curator"):
            _test(
                dict(
                    sub="not_test_user_id",
                    email="fake_user@email.com",
                ),
                additional_scope="write:collections",
            )

    def test__generate_s3_credentials__Not_Owner(self):
        collection = self.generate_collection(self.session, owner="not_test_user")
        response = self.app.post(
            f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials", headers=self.get_headers()
        )
        self.assertEqual(403, response.status_code, msg=response.data)

    def test__generate_s3_credentials__Not_Private(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name)
        response = self.app.post(
            f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials", headers=self.get_headers()
        )
        self.assertEqual(403, response.status_code)

    def get_headers(self):
        token = decode_token(get_auth_token(self.app)[8:].split(";")[0])
        return {"Authorization": f"Bearer {token['access_token']}"}


if __name__ == "__main__":
    unittest.main()
