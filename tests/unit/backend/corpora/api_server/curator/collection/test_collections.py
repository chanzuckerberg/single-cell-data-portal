import unittest
from unittest.mock import patch, Mock

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.lambdas.api.v1.authentication import decode_token
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token, make_token


class TestAuthToken(BaseAuthAPITest):
    @patch("backend.corpora.lambdas.api.v1.curation.collection.dataset.sts_client")
    def test__generate_s3_credentials__OK(self, sts_client: Mock):
        def _test(user_name: str, additional_scope: list = None):
            token_claims = dict(sub=user_name, email="fake_user@email.com")
            token = make_token(token_claims, additional_scope=additional_scope, token_duration=10)
            sts_client.assume_role_with_web_identity = Mock(
                return_value={
                    "Credentials": {
                        "AccessKeyId": "test_key",
                        "SecretAccessKey": "test_session_token",
                        "SessionToken": "test_session_token",
                    }
                }
            )
            collection = self.generate_collection(self.session)
            headers = {"Authorization": f"Bearer {token}"}

            response = self.app.post(
                f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials", headers=headers
            )
            self.assertEqual(200, response.status_code)
            self.assertEqual(response.json["Bucket"], "cellxgene-dataset-submissions-test")
            if additional_scope:
                self.assertEqual(response.json["UploadPath"], f"super_curator/{collection.id}/")
            else:
                self.assertEqual(response.json["UploadPath"], f"{user_name}/{collection.id}/")

        with self.subTest("collection owner"):
            _test(
                user_name="test_user_id",
            )

        with self.subTest("super curator"):
            _test(
                user_name="test_super_user_id",
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
