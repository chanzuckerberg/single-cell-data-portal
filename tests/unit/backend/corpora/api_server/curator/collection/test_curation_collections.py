import json
import unittest
from unittest.mock import patch, Mock

from backend.corpora.common.corpora_orm import CollectionVisibility, ProcessingStatus, DatasetArtifactFileType
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.corpora.api_server.mock_auth import make_token
from tests.unit.backend.fixtures.config import fake_s3_file


class TestAuthToken(BaseAuthAPITest):
    @patch("backend.corpora.lambdas.api.v1.curation.collections.collection_uuid.dataset.sts_client")
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
            self.assertEqual(response.json["UploadKeyPrefix"], f"{user_name}/{collection.id}/")

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
        tests = [
            (
                dict(
                    name="",
                    description="",
                    contact_name="",
                    contact_email="@email.com",
                    links=[{"link_type": "DOI", "link_url": "bad_doi"}],
                ),
                [
                    {"name": "contact_email", "reason": "Invalid format."},
                    {"name": "description", "reason": "Cannot be blank."},
                    {"name": "name", "reason": "Cannot be blank."},
                    {"name": "contact_name", "reason": "Cannot be blank."},
                    {"link_type": "DOI", "reason": "Invalid DOI"},
                ],
            ),
            (
                dict(
                    name="not blank",
                    description="description",
                    contact_name="some name",
                    contact_email="robot@email.com",
                    links=[
                        {"link_type": "DOI", "link_url": "doi:duplicated"},
                        {"link_type": "DOI", "link_url": "doi:duplicated"},
                    ],
                ),
                [{"link_type": "DOI", "reason": "Can only specify a single DOI"}],
            ),
        ]
        for body, expected_errors in tests:
            with self.subTest(body):
                response = self.app.post(
                    "/curation/v1/collections", headers=self.get_auth_headers(), data=json.dumps(body)
                )
                self.assertEqual(400, response.status_code)
                for error in expected_errors:
                    self.assertIn(error, response.json["detail"])


class TestGetCollections(BaseAuthAPITest):
    def setUp(self):
        super().setUp()
        self.test_collection = dict(
            name="collection", description="description", contact_name="john doe", contact_email="johndoe@email.com"
        )

    def test__get_collections_no_auth__OK(self):
        res_no_auth = self.app.get("/curation/v1/collections")
        self.assertEqual(200, res_no_auth.status_code)
        self.assertEqual(6, len(res_no_auth.json["collections"]))
        [self.assertEqual("PUBLIC", c["visibility"]) for c in res_no_auth.json["collections"]]

    def test__get_collections_with_auth__OK(self):
        res_auth = self.app.get("/curation/v1/collections", headers=self.get_auth_headers())
        self.assertEqual(200, res_auth.status_code)
        self.assertEqual(6, len(res_auth.json["collections"]))

    def test__get_collections_no_auth_visibility_private__OK(self):
        params = {"visibility": "PRIVATE"}
        res_private = self.app.get("/curation/v1/collections", query_string=params)
        self.assertEqual(401, res_private.status_code)

    def test__get_collections_no_auth_visibility_public__OK(self):
        params = {"visibility": "PUBLIC"}
        res_public = self.app.get("/curation/v1/collections", query_string=params)
        self.assertEqual(200, res_public.status_code)
        self.assertEqual(6, len(res_public.json["collections"]))

    def test__get_only_public_collections_with_auth__OK(self):
        params = {"visibility": "PUBLIC"}
        res = self.app.get("/curation/v1/collections", query_string=params, headers=self.get_auth_headers())
        self.assertEqual(200, res.status_code)
        self.assertEqual(6, len(res.json["collections"]))
        [self.assertEqual("PUBLIC", c["visibility"]) for c in res.json["collections"]]

    def test__get_only_private_collections_with_auth__OK(self):
        second_collection = self.generate_collection(self.session)
        for status in (ProcessingStatus.PENDING, ProcessingStatus.SUCCESS):
            self.generate_dataset(
                self.session,
                collection_id=second_collection.id,
                processing_status={"processing_status": status},
            ).id
        params = {"visibility": "PRIVATE"}
        res = self.app.get("/curation/v1/collections", query_string=params, headers=self.get_auth_headers())
        with self.subTest("Summary collection-level processing statuses are accurate"):
            for collection in res.json["collections"]:
                if collection["id"] == second_collection.id:
                    self.assertEqual(collection["processing_status"], "PENDING")
                else:
                    self.assertEqual(collection["processing_status"], "SUCCESS")
        self.assertEqual(200, res.status_code)
        self.assertEqual(2, len(res.json["collections"]))
        [self.assertEqual("PRIVATE", c["visibility"]) for c in res.json["collections"]]


class TestPutCollectionUUID(BaseAuthAPITest):
    def setUp(self):
        super().setUp()
        self.test_collection = dict(
            name="collection", description="description", contact_name="john doe", contact_email="johndoe@email.com"
        )
        self.generate_collection(
            self.session,
            id="test_curator_tag_collection_id",
            visibility=CollectionVisibility.PUBLIC.name,
            owner="test_user_id",
            name="test_collection_name",
            description="test_description",
            data_submission_policy_version="0",
            contact_name="Some Body",
            contact_email="somebody@chanzuckerberg.com",
        )
        self.generate_dataset(
            self.session,
            id="test_curator_tag",
            curator_tag="curator_tag",
            revision=0,
            name="test_dataset_name",
            schema_version="2.0.0",
            collection_id="test_curator_tag_collection_id",
            artifacts=[
                dict(
                    filename="test_filename",
                    filetype=DatasetArtifactFileType.H5AD.name,
                    user_submitted=True,
                    s3_uri=fake_s3_file,
                )
            ],
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
        headers = self.make_super_curator_header()
        response = self.app.put(
            f"/curation/v1/collections/{collection_uuid}", data=json.dumps(self.test_collection), headers=headers
        )
        self.assertEqual(200, response.status_code)


if __name__ == "__main__":
    unittest.main()
