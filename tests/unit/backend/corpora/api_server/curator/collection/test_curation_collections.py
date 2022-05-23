import json
import unittest
from unittest.mock import patch, Mock

from backend.corpora.common.corpora_orm import CollectionVisibility
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.corpora.api_server.mock_auth import make_token


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
            self.assertEqual(response.json["UploadKeyPrefix"], f"{user_name}/{collection.id}/")
            # if additional_scope:
            #     self.assertEqual(response.json["UploadKeyPrefix"], f"super|curator/{collection.id}/")
            # else:
            #     self.assertEqual(response.json["UploadKeyPrefix"], f"{user_name}/{collection.id}/")

        with self.subTest("collection owner"):
            _test(
                user_name="test_user_id",
            )

        # with self.subTest("super curator"):
        #     _test(
        #         user_name="test_super_user_id",
        #         additional_scope="write:collections",
        #     )

    def test__generate_s3_credentials__Not_Owner(self):
        collection = self.generate_collection(self.session, owner="not_test_user")
        response = self.app.post(
            f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials", headers=self.get_auth_headers()
        )
        self.assertEqual(403, response.status_code, msg=response.data)

    def test__generate_s3_credentials__Not_Private(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name)
        response = self.app.post(
            f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials", headers=self.get_auth_headers()
        )
        self.assertEqual(403, response.status_code)

    def test__generate_s3_credentials__No_Auth(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name)
        response = self.app.post(f"/curation/v1/collections/{collection.id}/datasets/s3-upload-credentials")
        self.assertEqual(401, response.status_code)


class TestPostCollection(BaseAuthAPITest):
    def setUp(self):
        super().setUp()
        self.test_collection = dict(
            name="collection", description="description", contact_name="john doe", contact_email="johndoe@email.com"
        )

    def test__create_collection__no_auth(self):
        response = self.app.post("/curation/v1/collections", data=json.dumps(self.test_collection))
        self.assertEqual(401, response.status_code)

    def test__create_collection__OK(self):
        response = self.app.post(
            "/curation/v1/collections", headers=self.get_auth_headers(), data=json.dumps(self.test_collection)
        )
        self.assertEqual(201, response.status_code)

    def test__create_collection__InvalidParameters(self):
        test_collection = dict(name="", description="", contact_name="", contact_email="@email.com")
        response = self.app.post(
            "/curation/v1/collections", headers=self.get_auth_headers(), data=json.dumps(test_collection)
        )
        self.assertEqual(400, response.status_code)


class TestPutCollectionUUID(BaseAuthAPITest):
    def setUp(self):
        super().setUp()
        self.test_collection = dict(
            name="collection", description="description", contact_name="john doe", contact_email="johndoe@email.com"
        )

    def test__update_collection__no_auth(self):
        collection_uuid = self.generate_collection(self.session).id
        response = self.app.put(f"/curation/v1/collections/{collection_uuid}", data=json.dumps(self.test_collection))
        self.assertEqual(401, response.status_code)

    def test__update_collection__OK(self):
        collection_uuid = self.generate_collection(self.session).id
        response = self.app.put(
            f"/curation/v1/collections/{collection_uuid}",
            data=json.dumps(self.test_collection),
            headers=self.get_auth_headers(),
        )
        self.assertEqual(200, response.status_code)

    def test__update_collection__Not_Owner(self):
        collection_uuid = self.generate_collection(self.session, owner="someone else").id
        response = self.app.put(
            f"/curation/v1/collections/{collection_uuid}",
            data=json.dumps(self.test_collection),
            headers=self.get_auth_headers(),
        )
        self.assertEqual(403, response.status_code)

    def test__update_collection__Super_Curator(self):
        collection_uuid = self.generate_collection(self.session).id
        headers = {"Authorization": "Bearer " + self.make_super_curator_token(), "Content-Type": "application/json"}
        response = self.app.put(
            f"/curation/v1/collections/{collection_uuid}", data=json.dumps(self.test_collection), headers=headers
        )
        self.assertEqual(200, response.status_code)


if __name__ == "__main__":
    unittest.main()
