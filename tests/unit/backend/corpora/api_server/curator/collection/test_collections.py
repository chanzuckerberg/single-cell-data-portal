import unittest
from unittest.mock import patch, Mock

from backend.corpora.common.corpora_orm import CollectionVisibility
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestAuthToken(BaseAuthAPITest):
    @patch("backend.corpora.lambdas.api.v1.authentication.assert_authorized_token")
    @patch("backend.corpora.lambdas.api.v1.curation.collection.dataset.sts_client")
    def test__generate_s3_credentials__OK(self, sts_client: Mock, assert_authorized_token: Mock):
        def _test(assert_authorized_token_return_value):
            assert_authorized_token.return_value = assert_authorized_token_return_value
            sts_client.assume_role_with_web_identity = Mock(
                return_value={
                    "access_key": "test_key",
                    "secret_access_key": "test_session_token",
                    "session_token": "test_session_token",
                }
            )
            collection = self.generate_collection(self.session)
            headers = {"Authorization": "Bearer fake_access_token", "IdToken": "test id token"}

            response = self.app.post(
                f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials", headers=headers
            )
            self.assertEqual(200, response.status_code)

        with self.subTest("collection owner"):
            _test({"sub": "test_user_id"})

        with self.subTest("super curator"):
            _test({"sub": "not_test_user_id", "scope": "write:collections"})

    @patch("backend.corpora.lambdas.api.v1.authentication.assert_authorized_token")
    def test__generate_s3_credentials__Not_Owner(self, assert_authorized_token: Mock):
        assert_authorized_token.return_value = {"sub": "test_user_id"}
        collection = self.generate_collection(self.session, owner="not_test_user")
        headers = {"Authorization": "Bearer fake_access_token"}
        response = self.app.post(
            f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials", headers=headers
        )
        self.assertEqual(403, response.status_code)

    @patch("backend.corpora.lambdas.api.v1.authentication.assert_authorized_token")
    def test__generate_s3_credentials__Not_Private(self, assert_authorized_token: Mock):
        assert_authorized_token.return_value = {"sub": "test_user_id"}
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name)
        headers = {"Authorization": "Bearer fake_access_token"}
        response = self.app.post(
            f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials", headers=headers
        )
        self.assertEqual(403, response.status_code)


if __name__ == "__main__":
    unittest.main()
