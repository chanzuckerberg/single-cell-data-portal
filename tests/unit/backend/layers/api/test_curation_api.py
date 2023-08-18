import copy
import json
import uuid
from collections import defaultdict
from dataclasses import asdict
from unittest.mock import Mock, patch

from backend.common.utils.api_key import generate
from backend.curation.api.v1.curation.collections.common import EntityColumns
from backend.layers.common.entities import (
    CollectionId,
    CollectionLinkType,
    CollectionMetadata,
    CollectionVersion,
    CollectionVisibility,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersionId,
    Link,
    OntologyTermId,
)
from backend.layers.thirdparty.crossref_provider import CrossrefDOINotFoundException
from tests.unit.backend.layers.api.test_portal_api import generate_mock_publisher_metadata
from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest
from tests.unit.backend.layers.common.base_test import DatasetArtifactUpdate, DatasetStatusUpdate

mock_config_attr = {
    "curator_role_arn": "test_role_arn",
    "submission_bucket": "cellxgene-dataset-submissions-test",
    "upload_max_file_size_gb": 1,
    "dataset_assets_base_url": "http://domain",
}


def mock_config_fn(name):
    return mock_config_attr[name]


class TestDeleteCollection(BaseAPIPortalTest):
    def _test(self, collection_id, header, expected_status, query_param_str=None):
        if header == "owner":
            headers = self.make_owner_header()
        elif header == "super":
            headers = self.make_super_curator_header()
        elif header == "cxg_admin":
            headers = self.make_cxg_admin_header()
        elif header == "not_owner":
            headers = self.make_not_owner_header()
        elif "noauth":
            headers = {}
        response = self.app.delete(
            f"/curation/v1/collections/{collection_id}{'?' + query_param_str if query_param_str else ''}",
            headers=headers,
        )
        self.assertEqual(expected_status, response.status_code)

    def test__delete_public_collection(self):
        tests = [
            ("not_owner", (403, 403)),
            ("noauth", (401, 401)),
            ("owner", (405, 405)),
            ("super", (405, 405)),
            ("cxg_admin", (405, 204)),
        ]
        query_param_strs = ("", "delete_published=true")
        public_collection_id = self.generate_published_collection().collection_id
        for auth, expected_responses in tests:
            for i, query_param_str in enumerate(query_param_strs):
                with self.subTest(auth):
                    self._test(public_collection_id, auth, expected_responses[i], query_param_str=query_param_str)

    def test__delete_revision_collection_by_collection_version_id(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 204), ("super", 204)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                revision_collection = self.generate_collection_revision()
                self._test(revision_collection.version_id, auth, expected_response)
                if expected_response == 204:
                    self._test(revision_collection.version_id, auth, 404)

    def test__delete_private_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 204), ("super", 204)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                private_collection_id = self.generate_unpublished_collection().collection_id
                self._test(private_collection_id, auth, expected_response)
                if expected_response == 204:
                    self._test(private_collection_id, auth, 404)

    def test__delete_private_collection_by_collection_version_id(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 403), ("super", 403), ("cxg_admin", 403)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                private_collection_version_id = self.generate_unpublished_collection().version_id
                self._test(private_collection_version_id, auth, expected_response)

    def test__delete_tombstone_collection(self):
        tests = [
            ("not_owner", (410, 410)),
            ("noauth", (401, 401)),
            ("owner", (410, 410)),
            ("super", (410, 410)),
            ("cxg_admin", (410, 410)),
        ]
        query_param_strs = ("", "delete_published=true")
        for auth, expected_responses in tests:
            for i, query_param_str in enumerate(query_param_strs):
                with self.subTest(auth):
                    collection = self.generate_published_collection()
                    self.business_logic.tombstone_collection(collection.collection_id)
                    self._test(collection.collection_id, auth, expected_responses[i], query_param_str=query_param_str)


class TestS3Credentials(BaseAPIPortalTest):
    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    @patch("backend.curation.api.v1.curation.collections.collection_id.s3_upload_credentials.sts_client")
    def test__generate_s3_credentials__OK(self, sts_client: Mock, mock_config: Mock):
        def _test(token, is_super_curator: bool = False):
            sts_client.assume_role_with_web_identity = Mock(
                return_value={
                    "Credentials": {
                        "AccessKeyId": "test_key",
                        "SecretAccessKey": "test_session_token",
                        "SessionToken": "test_session_token",
                    }
                }
            )
            headers = {"Authorization": f"Bearer {token}"}

            with self.subTest(f"{token}, unpublished collection"):
                unpublished_collection = self.generate_unpublished_collection()

                _id = unpublished_collection.collection_id
                response = self.app.get(f"/curation/v1/collections/{_id}/s3-upload-credentials", headers=headers)
                self.assertEqual(200, response.status_code)
                token_sub = self._mock_assert_authorized_token(token)["sub"]
                self.assertEqual(response.json["Bucket"], "cellxgene-dataset-submissions-test")
                if is_super_curator:
                    self.assertEqual(response.json["UploadKeyPrefix"], f"super/{_id}/")
                else:
                    self.assertEqual(response.json["UploadKeyPrefix"], f"{token_sub}/{_id}/")

            with self.subTest(f"{token}, unpublished revision"):
                collection_id = self.generate_published_collection().collection_id
                _id = self.generate_revision(collection_id).version_id
                response = self.app.get(f"/curation/v1/collections/{_id}/s3-upload-credentials", headers=headers)
                self.assertEqual(200, response.status_code)
                token_sub = self._mock_assert_authorized_token(token)["sub"]
                self.assertEqual(response.json["Bucket"], "cellxgene-dataset-submissions-test")
                if is_super_curator:
                    self.assertEqual(response.json["UploadKeyPrefix"], f"super/{_id}/")
                else:
                    self.assertEqual(response.json["UploadKeyPrefix"], f"{token_sub}/{_id}/")

        _test("owner")
        _test("super", is_super_curator=True)

    def test__generate_s3_credentials__Not_Owner(self):
        collection_id = self.generate_unpublished_collection(owner="not_test_user").collection_id
        response = self.app.get(
            f"/curation/v1/collections/{collection_id}/s3-upload-credentials", headers=self.make_owner_header()
        )
        self.assertEqual(403, response.status_code, msg=response.data)

    def test__generate_s3_credentials__Not_Private(self):
        collection_id = self.generate_published_collection().collection_id
        response = self.app.get(
            f"/curation/v1/collections/{collection_id}/s3-upload-credentials", headers=self.make_owner_header()
        )
        self.assertEqual(403, response.status_code)

    def test__generate_s3_credentials__No_Auth(self):
        collection_id = self.generate_unpublished_collection().collection_id
        response = self.app.get(f"/curation/v1/collections/{collection_id}/s3-upload-credentials")
        self.assertEqual(401, response.status_code)

    def test__generate_s3_credentials__by_collection_version_id(self):
        unpublished_collection = self.generate_unpublished_collection()
        version_id = unpublished_collection.version_id

        with self.subTest("collection owner"):
            token = "owner"
            headers = {"Authorization": f"Bearer {token}"}
            response = self.app.get(f"/curation/v1/collections/{version_id}/s3-upload-credentials", headers=headers)
            self.assertEqual(403, response.status_code)

        with self.subTest("super curator"):
            token = "super"
            headers = {"Authorization": f"Bearer {token}"}
            response = self.app.get(f"/curation/v1/collections/{version_id}/s3-upload-credentials", headers=headers)
            self.assertEqual(403, response.status_code)


class TestPostCollection(BaseAPIPortalTest):
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
            "/curation/v1/collections", headers=self.make_owner_header(), data=json.dumps(self.test_collection)
        )
        self.assertIn("collection_id", response.json.keys())
        self.assertEqual(201, response.status_code)

        # Check that the collection_id is the canonical collection ID
        collection_id = response.json["collection_id"]
        version = self.business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
        self.assertEqual(version.collection_id.id, collection_id)

        with self.subTest("Collection fields are correct"):
            for field, value in self.test_collection.items():
                self.assertEqual(value, getattr(version.metadata, field))

        with self.subTest("Curator name is set correctly"):
            self.assertEqual("First Last", version.curator_name)

    def test__create_collection__InvalidParameters(self):
        requests = [
            (
                dict(
                    name="",
                    description="",
                    contact_name="",
                    contact_email="@email.com",
                    doi="10.111/not_curie_reference_format",
                    consortia=["Not a valid consortia!"],
                ),
                [
                    {"name": "contact_email", "reason": "Invalid format."},
                    {"name": "description", "reason": "Cannot be blank."},
                    {"name": "name", "reason": "Cannot be blank."},
                    {"name": "contact_name", "reason": "Cannot be blank."},
                    {"name": "DOI", "reason": "DOI must be a CURIE reference."},
                    {"name": "consortia", "reason": "Invalid consortia."},
                ],
                6,
            ),
            (
                dict(
                    name="Nonprintable\x0acharacters\x0bare_NOT_allowed_in_name",
                    description="But\x0anonprintable\x0acharacters\x0bARE_allowed_in_description",
                    contact_name="And\x0ain_contact_name",
                    contact_email="somebody@email.com",
                    doi="10.111/not_curie_reference_format",
                ),
                [
                    {"name": "name", "reason": "Invalid characters detected."},
                    {"name": "DOI", "reason": "DOI must be a CURIE reference."},
                ],
                2,
            ),
        ]
        for body, expected_errors, num_expected_errors in requests:
            response = self.app.post(
                "/curation/v1/collections", headers=self.make_owner_header(), data=json.dumps(body)
            )
            self.assertEqual(400, response.status_code)
            for error in expected_errors:
                self.assertIn(error, response.json["invalid_parameters"])
            self.assertEqual(len(response.json["invalid_parameters"]), num_expected_errors)

    def test__create_collection__strip_string_fields(self):

        collection_metadata = dict(
            name="collection   ",
            description="   description",
            contact_name="  john doe  ",
            contact_email="  johndoe@email.com",
            consortia=["Consortia 1   "],
        )

        response = self.app.post(
            "/curation/v1/collections", headers=self.make_owner_header(), data=json.dumps(collection_metadata)
        )
        self.assertIn("collection_id", response.json.keys())
        self.assertEqual(201, response.status_code)

        # Check that the collection_id is the canonical collection ID
        collection_id = response.json["collection_id"]
        version = self.business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
        self.assertEqual(version.collection_id.id, collection_id)

        self.assertEqual(version.metadata.description, collection_metadata["description"].strip())
        self.assertEqual(version.metadata.name, collection_metadata["name"].strip())
        self.assertEqual(version.metadata.contact_name, collection_metadata["contact_name"].strip())
        self.assertEqual(version.metadata.contact_email, collection_metadata["contact_email"].strip())
        self.assertEqual(version.metadata.consortia, ["Consortia 1"])


class TestGetCollections(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.test_collection = dict(
            name="collection", description="description", contact_name="john doe", contact_email="johndoe@email.com"
        )
        self.expected_dataset_columns = EntityColumns.dataset_metadata_preview_cols + [
            "dataset_id",
            "dataset_version_id",
        ]
        self.expected_collection_columns = EntityColumns.collections_cols.copy()
        self.expected_collection_columns.remove("tombstone")
        self.expected_collection_columns.remove("owner")
        self.expected_collection_columns.append("processing_status")
        self.expected_collection_columns.append("collection_url")
        self.expected_collection_columns.append("doi")

    def _test_response(self, visibility=None, auth=False, status_code=200) -> dict:
        kwargs = {}
        if auth:
            kwargs["headers"] = self.make_owner_header()
        if visibility:
            kwargs["query_string"] = {"visibility": visibility}

        response = self.app.get("/curation/v1/collections", **kwargs)
        self.assertEqual(status_code, response.status_code)
        return response.json

    def test__get_collections_no_auth__OK(self):
        self.generate_unpublished_collection()
        self.generate_published_collection()
        response_body = self._test_response()
        self.assertEqual(1, len(response_body))
        self.assertEqual("PUBLIC", response_body[0]["visibility"])
        self.assertIsNone(response_body[0]["revising_in"])

    def test__get_collections_with_auth__OK_1(self):
        "The 'revising_in' attribute is None for unauthorized public collections"
        collection_id = self.generate_published_collection(owner="Someone Else").collection_id
        self.generate_revision(collection_id)
        resp = self._test_response(auth=True)
        self.assertEqual(1, len(resp))
        self.assertEqual(collection_id.id, resp[0]["collection_id"])
        self.assertIsNone(resp[0]["revising_in"])

    def test__get_collections_with_auth__OK_2(self):
        "The 'revising_in' attribute is None for collections which lack a revision"
        collection_id = self.generate_published_collection(owner="Someone Else").collection_id
        resp = self._test_response(auth=True)
        self.assertEqual(1, len(resp))
        self.assertEqual(collection_id.id, resp[0]["collection_id"])
        self.assertIsNone(resp[0]["revising_in"])

    def test__get_collections_with_auth__OK_3(self):
        "The 'revising_in' attribute is equal to the id of the revision Collection"
        collection_id = self.generate_published_collection().collection_id
        revision_id = self.generate_revision(collection_id).version_id
        resp = self._test_response(auth=True)
        self.assertEqual(1, len(resp))
        self.assertEqual(collection_id.id, resp[0]["collection_id"])
        self.assertEqual(revision_id.id, resp[0]["revising_in"])
        self.assertIsNone(resp[0]["revision_of"])

    def test__get_collections_with_auth__OK_4(self):
        "revision_of contains the id to their published counterpart"
        published_collection_id = self.generate_published_collection().collection_id
        collection_revision_id = self.generate_revision(published_collection_id).version_id
        resp = self._test_response(visibility="PRIVATE", auth=True)
        self.assertEqual(1, len(resp))
        self.assertEqual(collection_revision_id.id, resp[0]["collection_version_id"])
        self.assertEqual(published_collection_id.id, resp[0]["revision_of"])

    def test__get_collections_with_auth__OK_5(self):
        "revision_of contains None if the collection is public"
        published_collection_id = self.generate_published_collection().collection_id
        resp = self._test_response(visibility="PUBLIC", auth=True)
        self.assertEqual(1, len(resp))
        resp_collection = resp[0]
        self.assertEqual(published_collection_id.id, resp_collection["collection_id"])
        self.assertIsNone(resp_collection["revision_of"])

    def test__get_collections_with_auth__OK_6(self):
        "revision_of contains None if the collection is unpublished"
        self.generate_unpublished_collection()
        resp = self._test_response(visibility="PRIVATE", auth=True)
        self.assertEqual(1, len(resp))
        resp_collection = resp[0]
        self.assertIsNone(resp_collection["revision_of"])

    def test__get_collections_no_auth_visibility_private__403(self):
        self._test_response(visibility="PRIVATE", status_code=403)

    def test__get_collections_no_auth_visibility_public__OK(self):
        published_collection_id = self.generate_published_collection().collection_id
        self.generate_unpublished_collection()
        self.generate_revision(published_collection_id)
        resp = self._test_response(visibility="PUBLIC")
        self.assertEqual(1, len(resp))

    def test__get_a_curators_collections(self):
        curator_name = "John Smith"
        self.generate_published_collection(curator_name="Not Smith")
        self.generate_published_collection(curator_name=curator_name)
        self.generate_unpublished_collection(curator_name=curator_name)

        def _test(query_param, headers, expected_number_of_results):
            response = self.app.get("/curation/v1/collections", query_string=query_param, headers=headers)
            self.assertEqual(200, response.status_code)
            self.assertEqual(expected_number_of_results, len(response.json))
            for collection in response.json:
                self.assertEqual(curator_name, collection["curator_name"])

        visibilities = ["PUBLIC", "PRIVATE"]
        for visibility in visibilities:
            params = {"curator": curator_name, "visibility": visibility}
            with self.subTest(f"regular curators can't search by curator {visibility}"):
                response = self.app.get(
                    "/curation/v1/collections", query_string=params, headers=self.make_not_owner_header()
                )
                self.assertEqual(403, response.status_code)
            with self.subTest(f"users can't search by curator {visibility}"):
                response = self.app.get("/curation/v1/collections", query_string=params)
                self.assertEqual(403, response.status_code)

        with self.subTest("Searching for a curator that doesn't exists return 0 collections"):
            _test({"curator": "Not A Curator"}, self.make_super_curator_header(), 0)
        with self.subTest("super curators can search public collections by curator"):
            _test({"curator": curator_name}, self.make_super_curator_header(), 1)
        with self.subTest("super curators can search private collections by curator"):
            _test({"curator": curator_name, "visibility": "PRIVATE"}, self.make_super_curator_header(), 1)

    def test__get_only_private_collections_with_auth__OK(self):
        collection_version = self.generate_published_collection()
        self.generate_revision(collection_version.collection_id)
        self.generate_unpublished_collection()
        resp = self._test_response(visibility="PRIVATE", auth=True)
        self.assertEqual(2, len(resp))
        [self.assertEqual("PRIVATE", c["visibility"]) for c in resp]

    def test__collection_level_processing_status__INITIALIZED_as_PENDING(self):
        collection_version = self.generate_unpublished_collection(add_datasets=0)
        for status in (DatasetProcessingStatus.INITIALIZED, DatasetProcessingStatus.SUCCESS):
            self.generate_dataset(
                collection_version=collection_version,
                statuses=[DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=status)],
            )
        resp = self._test_response(visibility="PRIVATE", auth=True)
        self.assertEqual(1, len(resp))
        self.assertEqual("PENDING", resp[0]["processing_status"])

    def test__collection_level_processing_status__PENDING(self):
        collection_version = self.generate_unpublished_collection(add_datasets=0)
        for status in (DatasetProcessingStatus.PENDING, DatasetProcessingStatus.SUCCESS):
            self.generate_dataset(
                collection_version=collection_version,
                statuses=[DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=status)],
            )
        resp = self._test_response(visibility="PRIVATE", auth=True)
        self.assertEqual(1, len(resp))
        self.assertEqual("PENDING", resp[0]["processing_status"])

    def test__collection_level_processing_status__FAILURE(self):
        collection_version = self.generate_unpublished_collection(add_datasets=0)
        for status in (DatasetProcessingStatus.FAILURE, DatasetProcessingStatus.SUCCESS):
            self.generate_dataset(
                collection_version=collection_version,
                statuses=[DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=status)],
            )
        resp = self._test_response(visibility="PRIVATE", auth=True)
        self.assertEqual(1, len(resp))
        self.assertEqual("FAILURE", resp[0]["processing_status"])

    def test__collection_level_processing_status__SUCCESS(self):
        collection_version = self.generate_unpublished_collection(add_datasets=0)
        self.generate_dataset(
            collection_version=collection_version,
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.SUCCESS)
            ],
        )
        resp = self._test_response(visibility="PRIVATE", auth=True)
        self.assertEqual(1, len(resp))
        self.assertEqual("SUCCESS", resp[0]["processing_status"])

    def test__verify_expected_public_collection_fields(self):
        public_collection = self.generate_collection(
            visibility=CollectionVisibility.PUBLIC.name,
            links=[
                {
                    "link_name": "test_raw_data_link_name",
                    "link_type": "RAW_DATA",
                    "link_url": "http://test_raw_data_url.place",
                }
            ],
        )
        resp = self._test_response()
        self.assertEqual(1, len(resp))
        resp_collection = resp[0]
        self.check_fields(EntityColumns.link_cols, resp_collection["links"][0], "links")
        self.assertEqual(public_collection.datasets[0].dataset_id.id, resp_collection["datasets"][0]["dataset_id"])
        self.expected_collection_columns.remove("processing_status")
        self.check_fields(self.expected_dataset_columns, resp_collection["datasets"][0], "datasets")
        self.check_fields(self.expected_collection_columns, resp_collection, "collection")

    def test__verify_expected_private_collection_fields(self):
        private_collection = self.generate_collection(
            visibility=CollectionVisibility.PRIVATE.name,
            links=[
                {
                    "link_name": "test_raw_data_link_name",
                    "link_type": "RAW_DATA",
                    "link_url": "http://test_raw_data_url.place",
                }
            ],
            add_datasets=1,
        )

        body = self._test_response("PRIVATE", auth=True)
        self.assertEqual(1, len(body))
        resp_collection = body[0]

        self.check_fields(EntityColumns.link_cols, resp_collection["links"][0], "links")
        self.assertEqual(private_collection.datasets[0].dataset_id.id, resp_collection["datasets"][0]["dataset_id"])
        self.check_fields(self.expected_dataset_columns, resp_collection["datasets"][0], "datasets")
        self.check_fields(self.expected_collection_columns, resp_collection, "collection")

    def test__verify_expected_revision_collection_fields(self):
        published_version = self.generate_collection(
            visibility=CollectionVisibility.PUBLIC.name,
            links=[
                {
                    "link_name": "test_raw_data_link_name",
                    "link_type": "RAW_DATA",
                    "link_url": "http://test_raw_data_url.place",
                }
            ],
        )
        unpublished_version = self.generate_revision(published_version.collection_id)

        body = self._test_response("PRIVATE", auth=True)
        self.assertEqual(1, len(body))
        resp_collection = body[0]

        self.check_fields(EntityColumns.link_cols, resp_collection["links"][0], "links")
        self.assertEqual(
            unpublished_version.datasets[0].version_id.id, resp_collection["datasets"][0]["dataset_version_id"]
        )
        self.check_fields(self.expected_dataset_columns, resp_collection["datasets"][0], "datasets")
        self.check_fields(self.expected_collection_columns, resp_collection, "collection")

    def check_fields(self, fields: list, response: dict, entity: str):
        for key in fields:
            with self.subTest(f"{entity}:{key}"):
                self.assertIn(key, response.keys())
                response.pop(key)
        with self.subTest(f"No Extra fields in {entity}"):
            self.assertFalse(response)


class TestGetCollectionVersions(BaseAPIPortalTest):
    def confirm_timestamp_fields_present_then_remove(self, received_body: dict):
        # Confirm fields are present on Collection version body but ignore equality comparison for timestamps
        self.assertIn("created_at", received_body)
        self.assertIn("published_at", received_body)
        [self.assertIn("published_at", d) for d in received_body["dataset_versions"]]
        received_body.pop("created_at")
        received_body.pop("published_at")
        [d.pop("published_at") for d in received_body["dataset_versions"]]

    def test__get_collection_versions__200(self):
        # Create published collection with 2 published revisions and 1 unpublished revision
        published_collection = self.generate_published_collection()
        expected_version_ids = [published_collection.version_id.id]

        revision_collection = self.generate_revision(collection_id=published_collection.collection_id)
        expected_version_ids.append(revision_collection.version_id.id)
        self.business_logic.publish_collection_version(revision_collection.version_id)

        revision_collection_2 = self.generate_revision(collection_id=published_collection.collection_id)
        expected_version_ids.append(revision_collection_2.version_id.id)
        self.business_logic.publish_collection_version(revision_collection_2.version_id)
        expected_version_ids.reverse()

        with self.subTest("Published versions are returned in reverse chronological order with no revision open"):
            resp = self.app.get(f"/curation/v1/collections/{published_collection.collection_id.id}/versions")
            received_version_ids = [c_v["collection_version_id"] for c_v in resp.json]
            [self.confirm_timestamp_fields_present_then_remove(c_v) for c_v in resp.json]
            self.assertEqual(expected_version_ids, received_version_ids)

        self.generate_revision(collection_id=published_collection.collection_id)
        with self.subTest("Published versions are returned in reverse chronological order with a revision open"):
            resp = self.app.get(f"/curation/v1/collections/{published_collection.collection_id.id}/versions")
            received_version_ids = [c_v["collection_version_id"] for c_v in resp.json]
            [self.confirm_timestamp_fields_present_then_remove(c_v) for c_v in resp.json]
            self.assertEqual(expected_version_ids, received_version_ids)

    def test__get_collection_versions_tombstoned__410(self):
        published_collection = self.generate_published_collection()
        self.business_logic.tombstone_collection(published_collection.collection_id)
        with self.subTest("Returns 410 when a tombstoned canonical id is requested"):
            resp = self.app.get(f"/curation/v1/collections/{published_collection.collection_id.id}/versions")
            self.assertEqual(410, resp.status_code)

    def test__get_collection_versions_not_published_canonical__404(self):
        published_collection = self.generate_published_collection()
        revision_collection = self.generate_revision(collection_id=published_collection.collection_id)
        unpublished_collection = self.generate_unpublished_collection()
        with self.subTest("Returns 404 when nonexistent id is requested"):
            resp = self.app.get("/curation/v1/collections/01234567-89ab-cdef-0123-456789abcdef/versions")
            self.assertEqual(404, resp.status_code)

        with self.subTest("Returns 404 when revision id is requested"):
            resp = self.app.get(f"/curation/v1/collections/{revision_collection.version_id.id}/versions")
            self.assertEqual(404, resp.status_code)

        with self.subTest("Returns 404 when published version id is requested"):
            resp = self.app.get(f"/curation/v1/collections/{published_collection.version_id.id}/versions")
            self.assertEqual(404, resp.status_code)

        with self.subTest("Returns empty list when unpublished canonical id is requested"):
            resp = self.app.get(f"/curation/v1/collections/{unpublished_collection.collection_id.id}/versions")
            self.assertEqual(404, resp.status_code)

    def test__get_collection_versions_not_a_uuid__403(self):
        self.generate_published_collection()
        with self.subTest("Returns 403 when id is not uuid format"):
            resp = self.app.get("/curation/v1/collections/this_identifier_is_not_a_uuid/versions")
            self.assertEqual(403, resp.status_code)


class TestGetCollectionID(BaseAPIPortalTest):
    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    def test__get_published_collection_verify_body_is_reshaped_correctly__OK(self, mock_config: Mock):
        # Setup
        # test fixtures
        dataset_metadata = copy.deepcopy(self.sample_dataset_metadata)
        dataset_metadata.sex.append(OntologyTermId(label="test_sex2", ontology_term_id="test_obp"))
        links = [
            {
                "link_name": "test_raw_data_link_name",
                "link_type": "RAW_DATA",
                "link_url": "http://test_raw_data_url.place",
            },
            {"link_name": None, "link_type": "RAW_DATA", "link_url": "http://test_no_link_name_raw_data_url.place"},
            {
                "link_name": "test_protocol_link_name",
                "link_type": "PROTOCOL",
                "link_url": "http://test_protocol_url.place",
            },
            {"link_name": None, "link_type": "PROTOCOL", "link_url": "http://test_no_link_name_protocol_url.place"},
            {
                "link_name": "test_lab_website_link_name",
                "link_type": "LAB_WEBSITE",
                "link_url": "http://test_lab_website_url.place",
            },
            {
                "link_name": None,
                "link_type": "LAB_WEBSITE",
                "link_url": "http://test_no_link_name_lab_website_url.place",
            },
            {"link_name": "test_other_link_name", "link_type": "OTHER", "link_url": "http://test_other_url.place"},
            {"link_name": None, "link_type": "OTHER", "link_url": "http://test_no_link_name_other_url.place"},
            {
                "link_name": "test_data_source_link_name",
                "link_type": "DATA_SOURCE",
                "link_url": "http://test_data_source_url.place",
            },
            {
                "link_name": None,
                "link_type": "DATA_SOURCE",
                "link_url": "http://test_no_link_name_data_source_url.place",
            },
        ]
        collection_version = self.generate_collection(links=links, visibility="PRIVATE")
        self.generate_dataset(
            collection_version=collection_version,
            metadata=dataset_metadata,
            artifacts=[
                DatasetArtifactUpdate(type="h5ad", uri="http://test_filename/1234-5678-9/local.h5ad"),
                DatasetArtifactUpdate(type="rds", uri="http://test_filename/1234-5678-9/local.rds"),
                DatasetArtifactUpdate(type="cxg", uri="http://test_filename/1234-5678-9/local.cxg"),
                DatasetArtifactUpdate(type="raw_h5ad", uri="http://test_filename/1234-5678-9/raw.h5ad"),
            ],
        )
        self.business_logic.publish_collection_version(collection_version.version_id)
        collection_version = self.business_logic.get_collection_version(collection_version.version_id)
        dataset = collection_version.datasets[0]
        # expected results
        expect_dataset = asdict(collection_version.datasets[0].metadata)
        expect_dataset["title"] = expect_dataset.pop("name")
        expect_dataset.update(
            **{
                "explorer_url": f"/e/{dataset.dataset_id}.cxg/",
                "dataset_id": dataset.dataset_id.id,
                "dataset_version_id": dataset.version_id.id,
                "tombstone": False,
                "assets": [  # Filter out disallowed file types + properly construct url
                    {
                        "filesize": -1,
                        "filetype": "H5AD",
                        "url": f"http://domain/{dataset.version_id.id}.h5ad",
                    },
                    {
                        "filesize": -1,
                        "filetype": "RDS",
                        "url": f"http://domain/{dataset.version_id.id}.rds",
                    },
                ],
                "is_primary_data": [True, False],
                "x_approximate_distribution": "NORMAL",
            }
        )
        expected_body = asdict(collection_version.metadata)
        expected_body.update(
            **{
                "collection_id": collection_version.collection_id.id,
                "collection_url": f"https://frontend.corporanet.local:3000/collections/"
                f"{collection_version.collection_id}",
                "collection_version_id": collection_version.version_id.id,
                "datasets": [expect_dataset],
                "doi": None,
                "links": links,
                "publisher_metadata": None,
                "revision_of": None,
                "revising_in": None,
                "visibility": "PUBLIC",
                "curator_name": "Test User",
            }
        )

        # test
        res = self.app.get(f"/curation/v1/collections/{collection_version.collection_id}")
        self.assertEqual(200, res.status_code)
        res_body = res.json
        del res_body["created_at"]  # too finicky; ignore
        del res_body["revised_at"]  # too finicky; ignore
        del res_body["published_at"]  # too finicky; ignore
        for dataset in res_body["datasets"]:
            del dataset["revised_at"]  # too finicky; ignore
            del dataset["published_at"]  # too finicky; ignore
        self.maxDiff = None
        self.assertDictEqual(expected_body, res_body)  # Confirm dict has been packaged in list
        self.assertEqual(json.dumps(expected_body, sort_keys=True), json.dumps(res_body))

    def test__get_public_collection_verify_consortia_sorted__OK_1(self):
        collection_metadata = copy.deepcopy(self.sample_collection_metadata)
        collection_metadata.consortia = ["Consortia 3", "Consortia 1", "Consortia 2"]
        collection_version = self.generate_unpublished_collection(metadata=collection_metadata)
        self.assertEqual(collection_version.metadata.consortia, sorted(collection_metadata.consortia))

    def test__get_public_collection_verify_consortia_sorted__OK_2(self):
        collection_metadata = copy.deepcopy(self.sample_collection_metadata)
        collection_metadata.consortia = ["Consortia 3", "Consortia 1", "Consortia 2"]
        collection_version = self.generate_unpublished_collection(metadata=collection_metadata)

        res = self.app.get(f"/curation/v1/collections/{collection_version.collection_id}")
        self.assertEqual(res.json["consortia"], sorted(collection_metadata.consortia))

    def test__get_private_collection__OK(self):
        collection_version = self.generate_unpublished_collection()
        self._test_response(collection_version)

    def test__get_public_collection__OK(self):
        collection_version = self.generate_published_collection()
        self._test_response(collection_version)

    def test__get_unbpublished_revision__OK(self):
        collection_id = self.generate_published_collection().collection_id
        revision = self.generate_revision(collection_id)
        self._test_response(revision)

    def test_get_collection_dynamic_fields(self):
        def _test_responses(test_id, expected_id, auth_headers) -> dict:
            response = None
            for header in auth_headers:
                """Check that a request with different permissions levels returns same response."""
                resp = self.app.get(f"/curation/v1/collections/{test_id}", headers=header)
                self.assertEqual(200, resp.status_code)
                self.assertTrue(resp.json["collection_url"].endswith(expected_id.id))
                if not response:
                    response = resp
            return response.json

        privileged_access_headers = [self.make_owner_header(), self.make_super_curator_header()]
        restricted_access_headers = [self.make_not_owner_header(), self.make_not_auth_header()]
        all_headers = [*privileged_access_headers, *restricted_access_headers]
        unpublished = self.generate_unpublished_collection(add_datasets=1)
        with self.subTest("get unpublished version"):
            resp_collection = _test_responses(unpublished.collection_id, unpublished.collection_id, all_headers)
            self.assertEqual("PRIVATE", resp_collection["visibility"])
            self.assertIsNone(resp_collection.get("revising_in"))
            self.assertIsNone(resp_collection["revised_at"])
            self.assertIsNone(resp_collection["revision_of"])
            self.assertIsNone(resp_collection["published_at"])
            self.assertTrue(resp_collection["collection_url"].endswith(unpublished.collection_id.id))
            self.assertEqual(unpublished.collection_id.id, resp_collection["collection_id"])
            self.assertIsNone(resp_collection["datasets"][0]["revised_at"])
            self.assertEqual(unpublished.datasets[0].dataset_id.id, resp_collection["datasets"][0]["dataset_id"])
            self.assertIn(unpublished.datasets[0].dataset_id.id, resp_collection["datasets"][0]["explorer_url"])
            self.assertIn("processing_status", resp_collection["datasets"][0].keys())
            self.assertIn("processing_status", resp_collection.keys())

        published = self.generate_published_collection(add_datasets=1)
        with self.subTest("get published version"):
            resp_collection = _test_responses(published.collection_id, published.collection_id, all_headers)
            self.assertEqual("PUBLIC", resp_collection["visibility"])
            self.assertIsNone(resp_collection.get("revising_in"))
            self.assertIsNone(resp_collection["revised_at"])
            self.assertIsNone(resp_collection["revision_of"])
            self.assertIsNotNone(resp_collection["published_at"])
            self.assertTrue(resp_collection["collection_url"].endswith(published.collection_id.id))
            self.assertEqual(published.collection_id.id, resp_collection["collection_id"])
            self.assertIsNone(resp_collection["datasets"][0]["revised_at"])
            self.assertEqual(published.datasets[0].dataset_id.id, resp_collection["datasets"][0]["dataset_id"])
            self.assertIn(published.datasets[0].dataset_id.id, resp_collection["datasets"][0]["explorer_url"])
            self.assertNotIn("processing_status", resp_collection["datasets"][0].keys())
            self.assertNotIn("processing_status", resp_collection.keys())

        revision = self.generate_revision(published.collection_id)
        with self.subTest("get published with unpublished version and restricted access"):
            resp_collection = _test_responses(
                published.collection_id, published.collection_id, restricted_access_headers
            )
            self.assertEqual("PUBLIC", resp_collection["visibility"])
            self.assertIsNone(resp_collection.get("revising_in"))
            self.assertIsNone(resp_collection["revised_at"])
            self.assertIsNone(resp_collection["revision_of"])
            self.assertIsNotNone(resp_collection["published_at"])
            self.assertTrue(resp_collection["collection_url"].endswith(published.collection_id.id))
            self.assertEqual(published.collection_id.id, resp_collection["collection_id"])
            self.assertIsNone(resp_collection["datasets"][0]["revised_at"])
            self.assertEqual(published.datasets[0].dataset_id.id, resp_collection["datasets"][0]["dataset_id"])
            self.assertIn(published.datasets[0].dataset_id.id, resp_collection["datasets"][0]["explorer_url"])
            self.assertNotIn("processing_status", resp_collection["datasets"][0].keys())
            self.assertNotIn("processing_status", resp_collection.keys())

        with self.subTest("get published with unpublished version and privileged access"):
            resp_collection = _test_responses(
                published.collection_id, published.collection_id, privileged_access_headers
            )
            self.assertEqual("PUBLIC", resp_collection["visibility"])
            self.assertEqual(revision.version_id.id, resp_collection["revising_in"])
            self.assertIsNone(resp_collection["revised_at"])
            self.assertIsNone(resp_collection["revision_of"])
            self.assertIsNotNone(resp_collection["published_at"])
            self.assertTrue(resp_collection["collection_url"].endswith(published.collection_id.id))
            self.assertEqual(published.collection_id.id, resp_collection["collection_id"])
            self.assertIsNone(resp_collection["datasets"][0]["revised_at"])
            self.assertEqual(published.datasets[0].dataset_id.id, resp_collection["datasets"][0]["dataset_id"])
            self.assertIn(published.datasets[0].dataset_id.id, resp_collection["datasets"][0]["explorer_url"])
            self.assertNotIn("processing_status", resp_collection["datasets"][0].keys())
            self.assertNotIn("processing_status", resp_collection.keys())

        with self.subTest("get unpublished version (revision) with published version and read access"):
            resp_collection = _test_responses(revision.version_id, revision.version_id, all_headers)
            self.assertEqual("PRIVATE", resp_collection["visibility"])
            self.assertIsNone(resp_collection.get("revising_in"))
            self.assertEqual(revision.collection_id.id, resp_collection["revision_of"])
            self.assertIsNone(resp_collection["revised_at"])
            self.assertIsNotNone(resp_collection["published_at"])
            self.assertTrue(resp_collection["collection_url"].endswith(revision.version_id.id))
            self.assertEqual(revision.version_id.id, resp_collection["collection_version_id"])
            self.assertIsNone(resp_collection["datasets"][0]["revised_at"])
            self.assertEqual(revision.datasets[0].version_id.id, resp_collection["datasets"][0]["dataset_version_id"])
            self.assertIn(revision.datasets[0].version_id.id, resp_collection["datasets"][0]["explorer_url"])
            self.assertIn("processing_status", resp_collection["datasets"][0].keys())
            self.assertIn("processing_status", resp_collection.keys())

        revised_dataset = self.generate_dataset(
            collection_version=revision, replace_dataset_version_id=revision.datasets[0].version_id
        )
        with self.subTest("get unpublished version with replace dataset and published version"):
            resp_collection = _test_responses(revision.version_id, revision.version_id, all_headers)
            self.assertEqual("PRIVATE", resp_collection["visibility"])
            self.assertIsNone(resp_collection.get("revising_in"))
            self.assertEqual(revision.collection_id.id, resp_collection["revision_of"])
            self.assertIsNone(resp_collection["revised_at"])
            self.assertIsNotNone(resp_collection["published_at"])
            self.assertTrue(resp_collection["collection_url"].endswith(revision.version_id.id))
            self.assertEqual(revision.version_id.id, resp_collection["collection_version_id"])
            self.assertIsNone(resp_collection["datasets"][0]["revised_at"])
            self.assertEqual(revised_dataset.dataset_version_id, resp_collection["datasets"][0]["dataset_version_id"])
            self.assertIn(revised_dataset.dataset_version_id, resp_collection["datasets"][0]["explorer_url"])
            self.assertIn("processing_status", resp_collection["datasets"][0].keys())
            self.assertIn("processing_status", resp_collection.keys())

        self.business_logic.publish_collection_version(revision.version_id)
        with self.subTest("get updated published version"):
            resp_collection = _test_responses(published.collection_id, published.collection_id, all_headers)
            self.assertEqual("PUBLIC", resp_collection["visibility"])
            self.assertIsNone(resp_collection.get("revising_in"))
            self.assertIsNone(resp_collection["revision_of"])
            self.assertIsNotNone(resp_collection["revised_at"])
            self.assertIsNotNone(resp_collection["published_at"])
            self.assertTrue(resp_collection["collection_url"].endswith(published.collection_id.id))
            self.assertEqual(published.collection_id.id, resp_collection["collection_id"])
            self.assertIsNotNone(resp_collection["datasets"][0]["revised_at"])
            self.assertEqual(published.datasets[0].dataset_id.id, resp_collection["datasets"][0]["dataset_id"])
            self.assertIn(revision.datasets[0].dataset_id.id, resp_collection["datasets"][0]["explorer_url"])
            self.assertNotIn("processing_status", resp_collection["datasets"][0].keys())
            self.assertNotIn("processing_status", resp_collection.keys())

    def test__get_collection_with_dataset_failing_validation(self):
        collection_version = self.generate_collection(
            visibility=CollectionVisibility.PRIVATE.name,
        )
        dataset = self.generate_dataset(
            collection_version=collection_version,
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.FAILURE),
                DatasetStatusUpdate(status_key=DatasetStatusKey.VALIDATION, status=DatasetValidationStatus.INVALID),
            ],
            validation_message="test message",
        )
        res_json = self._test_response(collection_version)
        self.assertEqual("FAILURE", res_json["processing_status"])
        actual_dataset = res_json["datasets"][0]
        self.assertEqual(dataset.dataset_id, actual_dataset["dataset_id"])
        self.assertEqual("VALIDATION_FAILURE", actual_dataset["processing_status"])
        self.assertEqual("test message", actual_dataset["processing_status_detail"])

    def test__get_collection_with_dataset_failing_pipeline(self):
        collection = self.generate_collection(
            visibility=CollectionVisibility.PRIVATE.name,
        )
        dataset = self.generate_dataset(
            collection_version=collection,
            statuses=[
                DatasetStatusUpdate(status_key=DatasetStatusKey.PROCESSING, status=DatasetProcessingStatus.FAILURE)
            ],
        )
        res_json = self._test_response(collection)
        self.assertEqual("FAILURE", res_json["processing_status"])
        actual_dataset = res_json["datasets"][0]
        self.assertEqual(dataset.dataset_id, actual_dataset["dataset_id"])
        self.assertEqual("PIPELINE_FAILURE", actual_dataset["processing_status"])

    def test__get_nonexistent_collection__4xx(self):
        with self.subTest("Querying with historical collection version ID"):
            collection = self.generate_published_collection()
            revision = self.generate_revision(collection.collection_id)
            self.business_logic.publish_collection_version(revision.version_id)
            res = self.app.get(f"/curation/v1/collections/{collection.version_id}")
            self.assertEqual(403, res.status_code)
        with self.subTest("UUID input valid, but not found"):
            non_existent_id = str(uuid.uuid4())
            res = self.app.get(f"/curation/v1/collections/{non_existent_id}")
            self.assertEqual(404, res.status_code)
        with self.subTest("UUID input invalid"):
            non_existent_id = "123-example-fake-uuid"
            res = self.app.get(f"/curation/v1/collections/{non_existent_id}")
            self.assertEqual(403, res.status_code)

    def test__get_tombstoned_collection__410(self):
        collection_version = self.generate_published_collection()
        self.business_logic.tombstone_collection(collection_version.collection_id)
        self._test_response(collection_version, 410)

    def test_get_collection_with_no_datasets(self):
        collection_version = self.generate_unpublished_collection(add_datasets=0)
        self._test_response(collection_version)

    def test_get_collection_with_dataset_no_metadata(self):
        """
        GET collection should work when the collection has datasets with no metadata.
        This happens when the dataset did not complete ingestion yet.
        """
        collection_version = self.generate_unpublished_collection(add_datasets=0)
        self.business_logic.create_empty_dataset(collection_version.version_id)
        self._test_response(collection_version)

    def _test_response(self, collection_version: CollectionVersion, status_code=200, auth=True) -> dict:
        version_id = collection_version.version_id
        collection_id = collection_version.collection_id
        is_open_revision = (
            collection_version.published_at is None
            and collection_version.canonical_collection.originally_published_at is not None
        )
        headers = self.make_owner_header() if auth else self.make_not_owner_header()
        res_1 = self.app.get(f"/curation/v1/collections/{version_id}", headers=headers)
        if is_open_revision:
            self.assertEqual(status_code, res_1.status_code)
            if status_code == 200:
                self.assertEqual(collection_version.version_id.id, res_1.json["collection_id"])
        else:
            self.assertEqual(403, res_1.status_code)

        res_2 = self.app.get(f"/curation/v1/collections/{collection_id}", headers=headers)
        self.assertEqual(status_code, res_2.status_code)
        if status_code == 200:
            self.assertEqual(collection_version.collection_id.id, res_2.json["collection_id"])

        return res_1.json if is_open_revision else res_2.json

    def test__get_collection_with_x_approximate_distribution_none__OK(self):
        metadata = copy.deepcopy(self.sample_dataset_metadata)
        metadata.x_approximate_distribution = None
        dataset = self.generate_dataset(metadata=metadata, publish=True)
        res = self.app.get(f"/curation/v1/collections/{dataset.collection_id}", headers=self.make_owner_header())
        self.assertEqual(200, res.status_code)
        self.assertIsNone(res.json["datasets"][0]["x_approximate_distribution"])


class TestGetCollectionVersionID(BaseAPIPortalTest):
    def test_get_collection_version_ok(self):
        first_version = self.generate_published_collection()
        revision = self.generate_revision(first_version.collection_id)
        self.business_logic.publish_collection_version(revision.version_id)
        res = self.app.get(
            f"/curation/v1/collection_versions/{first_version.version_id}", headers=self.make_owner_header()
        )
        self.assertEqual(200, res.status_code)
        received_body = res.json
        expected_body = {
            "collection_id": f"{first_version.collection_id.id}",
            "collection_url": f"https://frontend.corporanet.local:3000/collections/{first_version.version_id.id}",
            "collection_version_id": f"{first_version.version_id.id}",
            "consortia": ["Consortia 1", "Consortia 2"],
            "contact_email": "john.doe@email.com",
            "contact_name": "john doe",
            "curator_name": "Jane Smith",
            "dataset_versions": [
                {
                    "assay": [{"label": "test_assay_label", "ontology_term_id": "test_assay_term_id"}],
                    "assets": [
                        {
                            "filesize": -1,
                            "filetype": "H5AD",
                            "url": f"None/{first_version.datasets[0].version_id.id}.h5ad",
                        },
                        {
                            "filesize": -1,
                            "filetype": "RDS",
                            "url": f"None/{first_version.datasets[0].version_id.id}.rds",
                        },
                    ],
                    "batch_condition": ["test_batch_1", "test_batch_2"],
                    "collection_id": f"{first_version.collection_id.id}",
                    "collection_version_id": f"{first_version.version_id.id}",
                    "cell_count": 10,
                    "cell_type": [{"label": "test_cell_type_label", "ontology_term_id": "test_cell_type_term_id"}],
                    "dataset_id": f"{first_version.datasets[0].dataset_id.id}",
                    "dataset_version_id": f"{first_version.datasets[0].version_id.id}",
                    "development_stage": [
                        {"label": "test_development_stage_label", "ontology_term_id": "test_development_stage_term_id"}
                    ],
                    "disease": [{"label": "test_disease_label", "ontology_term_id": "test_disease_term_id"}],
                    "donor_id": ["test_donor_1"],
                    "is_primary_data": [True, False],
                    "mean_genes_per_cell": 0.5,
                    "organism": [{"label": "test_organism_label", "ontology_term_id": "test_organism_term_id"}],
                    "schema_version": "3.0.0",
                    "self_reported_ethnicity": [
                        {
                            "label": "test_self_reported_ethnicity_label",
                            "ontology_term_id": "test_self_reported_ethnicity_term_id",
                        }
                    ],
                    "sex": [{"label": "test_sex_label", "ontology_term_id": "test_sex_term_id"}],
                    "suspension_type": ["test_suspension_type"],
                    "tissue": [{"label": "test_tissue_label", "ontology_term_id": "test_tissue_term_id"}],
                    "title": "test_dataset_name",
                    "tombstone": False,
                    "x_approximate_distribution": "NORMAL",
                }
            ],
            "description": "described",
            "doi": None,
            "links": [],
            "name": "test_collection",
            "publisher_metadata": None,
            "visibility": "PUBLIC",
        }
        # Confirm fields are present but ignore equality comparison for timestamps
        self.assertIn("created_at", received_body)
        received_body.pop("created_at")
        self.assertIn("published_at", received_body)
        received_body.pop("published_at")
        [self.assertIn("published_at", d) for d in received_body["dataset_versions"]]
        [d.pop("published_at") for d in received_body["dataset_versions"]]

        self.assertEqual(received_body, expected_body)

        res = self.app.get(f"/curation/v1/collection_versions/{revision.version_id}", headers=self.make_owner_header())
        self.assertEqual(200, res.status_code)
        self.assertEqual(res.json["collection_version_id"], revision.version_id.id)

    def test_get_collection_version_4xx(self):
        with self.subTest("Query endpoint with incorrect ID"):
            res = self.app.get(
                f"/curation/v1/collection_versions/{str(uuid.uuid4())}", headers=self.make_owner_header()
            )
            self.assertEqual(404, res.status_code)
        with self.subTest("Query endpoint with Canonical ID"):
            collection = self.generate_published_collection()
            res = self.app.get(
                f"/curation/v1/collection_versions/{collection.collection_id}", headers=self.make_owner_header()
            )
            self.assertEqual(404, res.status_code)
        with self.subTest("Query endpoint with non-UUID"):
            res = self.app.get("/curation/v1/collection_versions/bad-input-id", headers=self.make_owner_header())
            self.assertEqual(403, res.status_code)
        with self.subTest("Attempting to access tombstoned Collection via Collection version id returns 410 Gone"):
            collection = self.generate_published_collection()
            self.business_logic.tombstone_collection(collection.collection_id)
            res = self.app.get(
                f"/curation/v1/collection_versions/{collection.version_id}", headers=self.make_owner_header()
            )
            self.assertEqual(410, res.status_code)

        with self.subTest("Cannot access prior versions for a tombstoned Collection with multiple published versions"):
            collection = self.generate_published_collection()
            revision = self.generate_revision(collection.collection_id)
            self.business_logic.publish_collection_version(revision.version_id)
            self.business_logic.tombstone_collection(collection.collection_id)
            res = self.app.get(f"/curation/v1/collection_versions/{collection.version_id}")
            self.assertEqual(410, res.status_code)
            res = self.app.get(f"/curation/v1/collection_versions/{revision.version_id}")
            self.assertEqual(410, res.status_code)

        with self.subTest("Collection Version is unpublished collection"):
            collection = self.generate_unpublished_collection()
            res = self.app.get(
                f"/curation/v1/collection_versions/{collection.version_id}", headers=self.make_owner_header()
            )
            self.assertEqual(404, res.status_code)
        with self.subTest("Collection Version is unpublished revision"):
            first_version = self.generate_published_collection()
            revision = self.generate_revision(first_version.collection_id)
            res = self.app.get(
                f"/curation/v1/collection_versions/{revision.version_id}", headers=self.make_owner_header()
            )
            self.assertEqual(404, res.status_code)


class TestPatchCollectionID(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.test_collection = dict(
            name="collection", description="description", contact_name="john doe", contact_email="johndoe@email.com"
        )

    def test__update_collection__no_auth(self):
        collection_id = self.generate_collection(visibility="PRIVATE").collection_id
        response = self.app.patch(f"/curation/v1/collections/{collection_id}", data=json.dumps(self.test_collection))
        self.assertEqual(401, response.status_code)

    def test__update_collection__OK(self):
        collection_id = self.generate_collection(visibility="PRIVATE").collection_id
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            data=json.dumps(self.test_collection),
            headers=self.make_owner_header(),
        )
        self.assertEqual(200, response.status_code)

    def test__update_collection_by_version_id__403(self):
        version_id = self.generate_collection(visibility="PRIVATE").version_id
        response = self.app.patch(
            f"/curation/v1/collections/{version_id}",
            data=json.dumps(self.test_collection),
            headers=self.make_owner_header(),
        )
        self.assertEqual(403, response.status_code)

    def test__update_revision__OK(self):
        collection_id = self.generate_published_collection().collection_id
        revision_id = self.generate_revision(collection_id).version_id
        response = self.app.patch(
            f"/curation/v1/collections/{revision_id}",
            data=json.dumps(self.test_collection),
            headers=self.make_owner_header(),
        )
        self.assertEqual(200, response.status_code)

    def test__update_collection_partial_data__OK(self):
        links = [Link("name", "RAW_DATA", "http://test_link.place")]
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata())
        collection = self.generate_unpublished_collection(links=links)
        collection_id = collection.collection_id

        new_name = "A new name, and only a new name"
        metadata = {"name": new_name}
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            data=json.dumps(metadata),
            headers=self.make_owner_header(),
        )
        self.assertEqual(200, response.status_code)

        response = self.app.get(f"curation/v1/collections/{collection_id}")
        self.assertEqual(new_name, response.json["name"])
        self.assertEqual(collection.metadata.description, response.json["description"])
        self.assertEqual(collection.metadata.contact_name, response.json["contact_name"])
        self.assertEqual(collection.metadata.contact_email, response.json["contact_email"])
        self.assertEqual(
            [{"link_name": "name", "link_type": "RAW_DATA", "link_url": "http://test_link.place"}],
            response.json["links"],
        )
        self.assertEqual(collection.publisher_metadata, response.json["publisher_metadata"])

    def test__update_collection_strip_string_fields(self):
        links = [Link("   name ", "RAW_DATA", "http://test_link.place")]
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata())
        collection = self.generate_unpublished_collection(links=links)
        collection_id = collection.collection_id

        new_name = "A new name"
        new_description = "   description   "
        new_contact_name = "   jane dough"
        new_contact_email = "   janedough@email.com"
        new_consortia = ["Consortia 4   "]
        metadata = {
            "name": new_name,
            "description": new_description,
            "contact_name": new_contact_name,
            "contact_email": new_contact_email,
            "consortia": new_consortia,
        }
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            data=json.dumps(metadata),
            headers=self.make_owner_header(),
        )

        self.assertEqual(200, response.status_code)
        self.assertEqual(new_name.strip(), response.json["name"])
        self.assertEqual(new_description.strip(), response.json["description"])
        self.assertEqual(new_contact_name.strip(), response.json["contact_name"])
        self.assertEqual(new_contact_email.strip(), response.json["contact_email"])
        self.assertEqual(
            [{"link_name": "name", "link_type": "RAW_DATA", "link_url": "http://test_link.place"}],
            response.json["links"],
        )
        self.assertEqual(collection.publisher_metadata, response.json["publisher_metadata"])

        response = self.app.get(f"curation/v1/collections/{collection_id}")

        self.assertEqual(200, response.status_code)
        self.assertEqual(new_name.strip(), response.json["name"])
        self.assertEqual(new_description.strip(), response.json["description"])
        self.assertEqual(new_contact_name.strip(), response.json["contact_name"])
        self.assertEqual(new_contact_email.strip(), response.json["contact_email"])
        self.assertEqual(
            [{"link_name": "name", "link_type": "RAW_DATA", "link_url": "http://test_link.place"}],
            response.json["links"],
        )
        self.assertEqual(collection.publisher_metadata, response.json["publisher_metadata"])

    def test_update_collection_with_empty_required_fields(self):
        tests = [dict(description=""), dict(contact_name=""), dict(contact_email=""), dict(name="")]

        collection_id = self.generate_collection(visibility="PRIVATE").collection_id
        for test in tests:
            with self.subTest(test):
                response = self.app.patch(
                    f"/curation/v1/collections/{collection_id}",
                    data=json.dumps(test),
                    headers=self.make_owner_header(),
                )
                self.assertEqual(400, response.status_code)

    def test__update_collection__links__OK(self):
        links = [
            {"link_name": "name", "link_type": "RAW_DATA", "link_url": "http://test_link.place"},
        ]
        new_links = [
            {"link_name": "new link", "link_type": "RAW_DATA", "link_url": "http://brand_new_link.place"},
        ]

        links_configurations = (
            ("With links already in place; new links replace old", links, new_links, 200, new_links),
            ("With no links in place; new links get added", None, new_links, 200, new_links),
            ("With links in place, but empty request; no changes are made", links, None, 200, links),
            ("With links in place, empty array passed; BAD REQUEST 400", links, [], 400, links),
        )

        for test_title, initial_links, new_links, expected_status_code, expected_links in links_configurations:
            with self.subTest(test_title):
                collection_id = self.generate_collection(links=initial_links, visibility="PRIVATE").collection_id
                original_collection = self.app.get(f"curation/v1/collections/{collection_id}").json
                self.assertEqual(initial_links if initial_links else [], original_collection["links"])
                metadata = {"links": new_links} if new_links is not None else {}
                response = self.app.patch(
                    f"/curation/v1/collections/{collection_id}",
                    data=json.dumps(metadata),
                    headers=self.make_owner_header(),
                )
                self.assertEqual(expected_status_code, response.status_code)
                if expected_status_code == 200:
                    self.assertEqual(expected_links, response.json["links"])

    def test__update_collection__doi__OK(self):
        initial_doi = "10.2020"
        links = [
            {"link_name": "new doi", "link_type": "DOI", "link_url": "http://doi.org/10.2020"},
        ]
        new_doi = "10.1016"  # a real DOI (CURIE reference)
        collection_id = self.generate_collection(links=links, visibility="PRIVATE").collection_id
        original_collection = self.app.get(f"curation/v1/collections/{collection_id}").json
        self.assertEqual(initial_doi, original_collection["doi"])
        metadata = {"doi": new_doi}
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            json=metadata,
            headers=self.make_owner_header(),
        )
        self.assertEqual(200, response.status_code)
        self.assertEqual(new_doi, response.json["doi"])

    def test__update_collection__consortia__OK(self):
        initial_consortia = ["Consortia 1", "Consortia 2"]
        new_consortia = ["Consortia 3"]
        links = [
            {"link_name": "new doi", "link_type": "DOI", "link_url": "http://doi.org/10.2020"},
        ]
        collection_id = self.generate_collection(links=links, visibility="PRIVATE").collection_id
        original_collection = self.app.get(f"curation/v1/collections/{collection_id}").json
        self.assertEqual(initial_consortia, original_collection["consortia"])
        metadata = {"consortia": new_consortia}
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            json=metadata,
            headers=self.make_owner_header(),
        )
        self.assertEqual(200, response.status_code)
        self.assertEqual(new_consortia, response.json["consortia"])

    def test__remove_collection__consortia__OK(self):
        initial_consortia = ["Consortia 1", "Consortia 2"]
        new_consortia = []
        links = [
            {"link_name": "new doi", "link_type": "DOI", "link_url": "http://doi.org/10.2020"},
        ]
        collection_id = self.generate_collection(links=links, visibility="PRIVATE").collection_id
        original_collection = self.app.get(f"curation/v1/collections/{collection_id}").json
        self.assertEqual(initial_consortia, original_collection["consortia"])
        metadata = {"consortia": new_consortia}
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            json=metadata,
            headers=self.make_owner_header(),
        )
        self.assertEqual(200, response.status_code)
        self.assertEqual(new_consortia, response.json["consortia"])

    def test__update_public_collection_verify_fix_consortia_sort_order_OK(self):
        initial_consortia = ["Consortia 1", "Consortia 2"]
        new_consortia = ["Consortia 3", "Consortia 1", "Consortia 2"]
        links = [
            {"link_name": "new doi", "link_type": "DOI", "link_url": "http://doi.org/10.2020"},
        ]
        collection_id = self.generate_collection(links=links, visibility="PRIVATE").collection_id
        original_collection = self.app.get(f"curation/v1/collections/{collection_id}").json
        self.assertEqual(initial_consortia, original_collection["consortia"])
        metadata = {"consortia": new_consortia}
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            json=metadata,
            headers=self.make_owner_header(),
        )
        self.assertEqual(200, response.status_code)
        self.assertEqual(sorted(new_consortia), response.json["consortia"])

    def test__update_collection__doi_is_not_CURIE_reference__BAD_REQUEST(self):
        links = [
            {"link_name": "doi", "link_type": "DOI", "link_url": "http://doi.doi/10.1011/something"},
        ]
        collection = self.generate_collection(links=links, visibility="PRIVATE")
        collection_id = collection.collection_id
        original_collection = self.app.get(f"curation/v1/collections/{collection_id}").json

        metadata = {"doi": "https://doi.org/10.1016"}
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            json=metadata,
            headers=self.make_owner_header(),
        )
        self.assertEqual(400, response.status_code)
        original_collection_unchanged = self.app.get(f"curation/v1/collections/{collection_id}").json
        self.assertEqual(original_collection["doi"], original_collection_unchanged["doi"])

    def test__update_collection__links_None_does_not_remove_publisher_metadata(self):
        links = [
            {"link_name": "doi", "link_type": "DOI", "link_url": "http://doi.doi/10.1011/something"},
        ]

        mock_publisher_metadata = generate_mock_publisher_metadata()
        self.crossref_provider.fetch_metadata = Mock(return_value=mock_publisher_metadata)

        collection = self.generate_collection(links=links, visibility="PRIVATE")
        collection_id = collection.collection_id
        original_collection = self.app.get(f"curation/v1/collections/{collection_id}").json
        self.assertIsNotNone(original_collection["publisher_metadata"])

        metadata = {"name": "new collection title"}
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            json=metadata,
            headers=self.make_owner_header(),
        )
        self.assertEqual(200, response.status_code)
        updated_collection = self.app.get(f"curation/v1/collections/{collection_id}").json
        self.assertIsNotNone(updated_collection["publisher_metadata"])
        self.assertEqual(updated_collection["links"], original_collection["links"])

    def test__update_collection__doi_does_not_exist__BAD_REQUEST(self):
        links = [
            {"link_name": "name", "link_type": "RAW_DATA", "link_url": "http://test_link.place"},
            {"link_name": "doi", "link_type": "DOI", "link_url": "http://doi.doi/10.1011/something"},
        ]
        new_links = [
            {"link_name": "new link", "link_type": "RAW_DATA", "link_url": "http://brand_new_link.place"},
        ]
        mock_publisher_metadata = generate_mock_publisher_metadata()
        self.crossref_provider.fetch_metadata = Mock(return_value=mock_publisher_metadata)
        collection = self.generate_collection(links=links, visibility="PRIVATE")
        self.assertIsNotNone(collection.publisher_metadata)
        collection_id = collection.collection_id
        original_collection = self.app.get(f"curation/v1/collections/{collection_id}").json

        self.crossref_provider.fetch_metadata = Mock(side_effect=CrossrefDOINotFoundException())

        # Only compare to first item in links list because "DOI" type gets removed from Curator API response
        self.assertEqual(links[:1], original_collection["links"])
        metadata = {"links": new_links, "doi": "10.1016/bad_doi"}
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            data=json.dumps(metadata),
            headers=self.make_owner_header(),
        )
        self.assertEqual(400, response.status_code)
        original_collection_unchanged = self.app.get(f"curation/v1/collections/{collection_id}").json

        # Only compare to first item in links list because "DOI" type gets removed from Curator API response
        self.assertEqual(links[:1], original_collection_unchanged["links"])
        self.assertEqual(mock_publisher_metadata, original_collection_unchanged["publisher_metadata"])

    def test__update_collection__Not_Owner(self):
        collection_id = self.generate_unpublished_collection(owner="someone else").collection_id
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            data=json.dumps(self.test_collection),
            headers=self.make_owner_header(),
        )
        self.assertEqual(403, response.status_code)

    def test__update_collection__Super_Curator(self):
        collection_id = self.generate_unpublished_collection().collection_id
        headers = self.make_super_curator_header()
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}", data=json.dumps(self.test_collection), headers=headers
        )
        self.assertEqual(200, response.status_code)

    def test__update_public_collection_owner__405(self):
        collection_id = self.generate_published_collection().collection_id
        headers = self.make_super_curator_header()
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}", data=json.dumps(self.test_collection), headers=headers
        )
        self.assertEqual(405, response.status_code)
        self.assertEqual(
            "Directly editing a public Collection is not allowed; you must create a revision.",
            response.json["detail"],
        )


class TestDeleteDataset(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.auth_credentials = [
            (self.make_super_curator_header, "super", 202),
            (self.make_owner_header, "owner", 202),
            (None, "none", 401),
            (self.make_not_owner_header, "not_owner", 403),
        ]

    def _delete(self, auth, collection_id, dataset_id, query_param_str=None):
        """
        Helper method to call the delete endpoint
        """
        test_url = f"/curation/v1/collections/{collection_id}/datasets/{dataset_id}{'?' + query_param_str if query_param_str else ''}"
        headers = auth() if callable(auth) else auth
        return self.app.delete(test_url, headers=headers)

    def test__delete_dataset_by_version_id(self):
        """
        Calling DELETE /collections/:collection_id/datasets/:dataset_id should return a 403 in all cases with auth
        credentials, and a 401 without auth credentials
        """
        for auth, auth_description, _ in self.auth_credentials:
            with self.subTest(f"{auth_description}"):
                dataset = self.generate_dataset(
                    statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)],
                    publish=False,
                )
                response = self._delete(auth, dataset.collection_version_id, dataset.dataset_version_id)
                if auth:
                    self.assertEqual(403, response.status_code)
                else:
                    self.assertEqual(401, response.status_code)

    def test__delete_dataset_by_canonical_id(self):
        """
        Calling DELETE /collections/:collection_id/datasets/:dataset_id should work according to the
        auth token passed and when using canonical ids. In this case, the unpublished collection
        version will be looked up and used for deletion.
        """
        for auth, auth_description, expected_status_code in self.auth_credentials:
            with self.subTest(f"{auth_description} {expected_status_code}"):
                dataset = self.generate_dataset(
                    statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)],
                    publish=False,
                )
                response = self._delete(auth, dataset.collection_id, dataset.dataset_id)
                self.assertEqual(expected_status_code, response.status_code)

    def test__delete_published_dataset__405(self):
        """
        A Dataset that has been published cannot be deleted via the API
        """
        for auth_func in (self.make_super_curator_header, self.make_owner_header):
            with self.subTest("Cannot delete published Dataset in a revision"):
                collection = self.generate_published_collection()
                revision = self.generate_revision(collection.collection_id)
                dataset_id, dataset_version_id = revision.datasets[0].dataset_id, revision.datasets[0].version_id
                response = self._delete(auth_func, revision.collection_id, dataset_id)
                self.assertEqual(405, response.status_code)
            with self.subTest("Cannot delete published Dataset in a revision even after it has been updated"):
                self.generate_dataset(collection_version=revision, replace_dataset_version_id=dataset_version_id)
                response = self._delete(auth_func, revision.collection_id, dataset_id)
                self.assertEqual(405, response.status_code)

    def test__delete_published_dataset_cxg_admin(self):
        """
        cxg_admin role can delete published Datasets during revisions using query_param flag 'delete_published=true'
        """
        with self.subTest("Can delete published Dataset during a revision"):
            collection = self.generate_published_collection()
            revision = self.generate_revision(collection.collection_id)
            dataset_id = revision.datasets[0].dataset_id
            response = self._delete(
                self.make_cxg_admin_header, revision.version_id, dataset_id, query_param_str="delete_published=true"
            )
            self.assertEqual(202, response.status_code)
        with self.subTest("Cannot delete published Dataset without query param"):
            collection = self.generate_published_collection()
            revision = self.generate_revision(collection.collection_id)
            dataset_id = revision.datasets[0].dataset_id
            response = self._delete(self.make_cxg_admin_header, revision.version_id, dataset_id)
            self.assertEqual(405, response.status_code)
        with self.subTest("Cannot delete published Dataset from public Collection"):
            collection = self.generate_published_collection()
            revision = self.generate_revision(collection.collection_id)
            dataset_id = revision.datasets[0].dataset_id
            response = self._delete(
                self.make_cxg_admin_header, revision.collection_id, dataset_id, query_param_str="delete_published=true"
            )
            self.assertEqual(405, response.status_code)


class TestGetDatasets(BaseAPIPortalTest):
    def test_get_dataset_in_a_collection(self):
        dataset = self.generate_dataset(name="test")

        with self.subTest("by canonical dataset ID"):
            test_url = f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_id}"

            response = self.app.get(test_url)
            self.assertEqual(200, response.status_code)
            self.assertEqual(dataset.dataset_id, response.json["dataset_id"])

        with self.subTest("by dataset version ID"):
            test_url = f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}"

            response = self.app.get(test_url)
            self.assertEqual(403, response.status_code)

        collection = self.generate_published_collection()
        collection_id = collection.collection_id
        dataset_id = collection.datasets[0].dataset_id.id
        revision_id = self.generate_revision(collection_id).version_id
        self.business_logic.publish_collection_version(revision_id)
        unbpublished_revision_id = self.generate_revision(collection_id).version_id
        with self.subTest("by historical collection version ID"):
            test_url = f"/curation/v1/collections/{revision_id}/datasets/{dataset_id}"
            response = self.app.get(test_url)
            self.assertEqual(403, response.status_code)
        with self.subTest("by unpublished, active revision ID"):
            test_url = f"/curation/v1/collections/{unbpublished_revision_id}/datasets/{dataset_id}"
            response = self.app.get(test_url)
            self.assertEqual(200, response.status_code)
            self.assertEqual(dataset_id, response.json["dataset_id"])

    def test_get_dataset_in_a_tombstoned_collection_410(self):
        collection = self.generate_published_collection(add_datasets=2)
        dataset = collection.datasets[0]
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.dataset_id}"
        response = self.app.get(test_url)
        self.assertEqual(200, response.status_code)
        self.business_logic.tombstone_collection(collection.collection_id)
        response = self.app.get(test_url)
        self.assertEqual(410, response.status_code)

    def test_get_tombstoned_dataset_in_a_collection_410(self):
        collection = self.generate_published_collection(add_datasets=2)
        self.assertEqual(2, len(collection.datasets))
        dataset = collection.datasets[0]
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.dataset_id}"
        response = self.app.get(test_url)
        self.assertEqual(200, response.status_code)
        revision = self.generate_revision(collection.collection_id)
        self.business_logic.remove_dataset_version(revision.version_id, dataset.version_id, delete_published=True)
        self.business_logic.publish_collection_version(revision.version_id)
        new_published_version = self.database_provider.get_collection_version(revision.version_id)
        self.assertEqual(1, len(new_published_version.datasets))
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.dataset_id}"
        response = self.app.get(test_url)
        self.assertEqual(410, response.status_code)

    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    def test_get_dataset_shape(self, mock_config: Mock):
        # retrieve a private dataset
        private_dataset = self.generate_dataset(name="test")
        test_url = f"/curation/v1/collections/{private_dataset.collection_id}/datasets/{private_dataset.dataset_id}"
        response = self.app.get(test_url)
        body = response.json
        self.assertEqual("test", body["title"])
        expected_assets = [  # Filter out disallowed file types + properly construct url
            {
                "filesize": -1,
                "filetype": "H5AD",
                "url": f"http://domain/{private_dataset.dataset_version_id}.h5ad",
            },
            {
                "filesize": -1,
                "filetype": "RDS",
                "url": f"http://domain/{private_dataset.dataset_version_id}.rds",
            },
        ]
        self.assertEqual(expected_assets, body["assets"])

        # retrieve a public dataset
        public_dataset = self.generate_dataset(name="test", publish=True)
        test_url = f"/curation/v1/collections/{public_dataset.collection_id}/datasets/{public_dataset.dataset_id}"
        response = self.app.get(test_url)
        body = response.json
        self.assertEqual("test", body["title"])
        expected_assets = [  # Filter out disallowed file types + properly construct url
            {
                "filesize": -1,
                "filetype": "H5AD",
                "url": f"http://domain/{public_dataset.dataset_version_id}.h5ad",
            },
            {
                "filesize": -1,
                "filetype": "RDS",
                "url": f"http://domain/{public_dataset.dataset_version_id}.rds",
            },
        ]
        self.assertEqual(expected_assets, body["assets"])

        # retrieve a revised dataset using version_id
        collection_id = self.generate_published_collection(add_datasets=2).canonical_collection.id
        version = self.generate_revision(collection_id)
        dataset_version = self.generate_dataset(
            collection_version=version, replace_dataset_version_id=version.datasets[0].version_id
        )
        test_url = f"/curation/v1/collections/{version.version_id}/datasets/{dataset_version.dataset_id}"
        response = self.app.get(test_url)
        body = response.json
        expected_assets = [  # Filter out disallowed file types + properly construct url
            {
                "filesize": -1,
                "filetype": "H5AD",
                "url": f"http://domain/{dataset_version.dataset_version_id}.h5ad",
            },
            {
                "filesize": -1,
                "filetype": "RDS",
                "url": f"http://domain/{dataset_version.dataset_version_id}.rds",
            },
        ]
        self.assertEqual(expected_assets, body["assets"])
        self.assertEqual("test_dataset_name", body["title"])

        # retrieve an unrevised dataset in a revision Collection
        unreplaced_dataset = version.datasets[1]
        test_url = f"/curation/v1/collections/{version.version_id}/datasets/{unreplaced_dataset.dataset_id}"
        response = self.app.get(test_url)
        body = response.json
        self.assertEqual("test_dataset_name", body["title"])
        expected_assets = [  # Filter out disallowed file types + properly construct url
            {
                "filesize": -1,
                "filetype": "H5AD",
                "url": f"http://domain/{unreplaced_dataset.version_id}.h5ad",
            },
            {
                "filesize": -1,
                "filetype": "RDS",
                "url": f"http://domain/{unreplaced_dataset.version_id}.rds",
            },
        ]
        self.assertEqual(expected_assets, body["assets"])

        # retrieve a newly added dataset in a revision Collection
        new_dataset = self.generate_dataset(collection_version=version)
        test_url = f"/curation/v1/collections/{version.version_id}/datasets/{new_dataset.dataset_id}"
        response = self.app.get(test_url)
        body = response.json
        self.assertEqual("test_dataset_name", body["title"])
        expected_assets = [  # Filter out disallowed file types + properly construct url
            {
                "filesize": -1,
                "filetype": "H5AD",
                "url": f"http://domain/{new_dataset.dataset_version_id}.h5ad",
            },
            {
                "filesize": -1,
                "filetype": "RDS",
                "url": f"http://domain/{new_dataset.dataset_version_id}.rds",
            },
        ]
        self.assertEqual(expected_assets, body["assets"])

    def test_get_dataset_is_primary_data_shape(self):
        tests = [
            ("PRIMARY", [True]),
            ("SECONDARY", [False]),
            ("BOTH", [True, False]),
        ]
        for is_primary_data, result in tests:
            with self.subTest(f"{is_primary_data}=={result}"):
                metadata = self.sample_dataset_metadata
                metadata.is_primary_data = is_primary_data
                dataset = self.generate_dataset(metadata=metadata)
                test_url = f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_id}"
                response = self.app.get(test_url)
                self.assertEqual(result, response.json["is_primary_data"])

    def test_get_nonexistent_dataset_4xx(self):
        collection = self.generate_unpublished_collection(add_datasets=1)
        with self.subTest("Querying with version ID instead of canonical"):
            dataset = collection.datasets[0]
            test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.version_id}"
            response = self.app.get(test_url)
            self.assertEqual(403, response.status_code)
        with self.subTest("UUID input invalid"):
            non_existent_dataset_id = "123-example-fake-uuid"
            test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{non_existent_dataset_id}"
            response = self.app.get(test_url)
            self.assertEqual(403, response.status_code)
        with self.subTest("UUID input valid, but not found"):
            non_existent_dataset_id = str(uuid.uuid4())
            test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{non_existent_dataset_id}"
            response = self.app.get(test_url)
            self.assertEqual(404, response.status_code)

    def test_get_datasets_nonexistent_collection_4xx(self):
        with self.subTest("UUID input valid, but not found"):
            non_existent_id = str(uuid.uuid4())
            test_url = f"/curation/v1/collections/{non_existent_id}/datasets/{non_existent_id}"
            headers = self.make_owner_header()
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(404, response.status_code)
        with self.subTest("UUID input invalid"):
            non_existent_id = "123-example-fake-uuid"
            test_url = f"/curation/v1/collections/{non_existent_id}/datasets/{non_existent_id}"
            headers = self.make_owner_header()
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(403, response.status_code)

    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    def test_get_dataset_no_assets(self, mock_config: Mock):
        private_dataset = self.generate_dataset(artifacts=[])
        test_url = f"/curation/v1/collections/{private_dataset.collection_id}/datasets/{private_dataset.dataset_id}"
        response = self.app.get(test_url)
        body = response.json
        self.assertEqual([], body["assets"])

    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    def test_get_all_datasets_200(self, mock_config: Mock):
        published_collection_1 = self.generate_published_collection(
            add_datasets=2,
            metadata=CollectionMetadata(
                "test_collection_1",
                "described",
                "john doe",
                "john.doe@email.com",
                [Link(name="doi link", type=CollectionLinkType.DOI.name, uri="http://doi.org/12.3456/j.celrep")],
                ["Consortia 1", "Consortia 2"],
            ),
        )
        published_collection_2 = self.generate_published_collection(
            owner="other owner",
            curator_name="other curator",
            add_datasets=1,
            metadata=CollectionMetadata(
                "test_collection_2",
                "described",
                "john doe",
                "john.doe@email.com",
                [Link(name="doi link", type=CollectionLinkType.DOI.name, uri="http://doi.org/78.91011/j.celrep")],
                ["Consortia 1", "Consortia 2"],
            ),
        )
        self.generate_unpublished_collection(add_datasets=4)
        revision = self.generate_revision(published_collection_1.collection_id)

        with self.subTest("With super curator credentials"):
            headers = self.make_super_curator_header()
            super_curator_response = self.app.get("/curation/v1/datasets", headers=headers)
            self.assertEqual(3, len(super_curator_response.json))

        response = self.app.get("/curation/v1/datasets")

        with self.subTest("With no credentials"):
            self.assertEqual(3, len(response.json))

        with self.subTest("Contains collection_name and collection_doi"):
            collection_names = {published_collection_1.metadata.name, published_collection_2.metadata.name}
            expected_collection_dois = {"12.3456/j.celrep", "78.91011/j.celrep"}

            received_collection_names = set()
            received_collection_dois = set()
            for dataset in response.json:
                received_collection_names.add(dataset["collection_name"])
                received_collection_dois.add(dataset["collection_doi"])

            self.assertEqual(collection_names, received_collection_names)
            self.assertEqual(expected_collection_dois, received_collection_dois)

        self.generate_dataset(
            collection_version=revision,
            replace_dataset_version_id=revision.datasets[0].version_id,
            publish=True,
        )

        with self.subTest("Only public datasets are returned, in reverse-chronological order"):
            # Endpoint uses secondary sort by dataset_id for consistency with Datasets in same Collection
            sorted_dataset_ids = [published_collection_2.datasets[0].dataset_id.id] + sorted(
                [d.dataset_id.id for d in published_collection_1.datasets], reverse=True
            )
            received_dataset_ids = []
            for dataset in response.json:
                received_dataset_ids.append(dataset["dataset_id"])

            self.assertEqual(sorted_dataset_ids, received_dataset_ids)

        with self.subTest("The 'revised_at' field is populated for revised Datasets only"):
            response = self.app.get(
                "/curation/v1/datasets"
            )  # Refresh response because 'revised_at' is calculated field
            for dataset in response.json:
                if dataset["dataset_id"] == published_collection_1.datasets[0].dataset_id.id:
                    self.assertIsNotNone(dataset["revised_at"])
                else:
                    self.assertIsNone(dataset["revised_at"])

        with self.subTest("Response Dataset objects contain index-specific attributes"):
            index_specific_attributes = ("collection_doi", "collection_id", "collection_name")
            dataset = response.json[0]
            for attribute in index_specific_attributes:
                self.assertIn(attribute, dataset)

        with self.subTest("Response Dataset objects contain only expected assets"):
            dataset = response.json[0]
            expected_assets = [  # Filter out disallowed file types + properly construct url
                {
                    "filesize": -1,
                    "filetype": "H5AD",
                    "url": f"http://domain/{dataset['dataset_version_id']}.h5ad",
                },
                {
                    "filesize": -1,
                    "filetype": "RDS",
                    "url": f"http://domain/{dataset['dataset_version_id']}.rds",
                },
            ]
            self.assertEqual(expected_assets, dataset["assets"])

    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    def test_get_datasets_by_schema_200(self, mock_config: Mock):
        published_collection_1 = self.generate_published_collection(
            add_datasets=2,
            metadata=CollectionMetadata(
                "test_collection_1",
                "described",
                "john doe",
                "john.doe@email.com",
                [Link(name="doi link", type=CollectionLinkType.DOI.name, uri="http://doi.org/12.3456/j.celrep")],
                ["Consortia 1", "Consortia 2"],
            ),
        )
        published_collection_2 = self.generate_published_collection(
            owner="other owner",
            curator_name="other curator",
            add_datasets=1,
            metadata=CollectionMetadata(
                "test_collection_2",
                "described",
                "john doe",
                "john.doe@email.com",
                [Link(name="doi link", type=CollectionLinkType.DOI.name, uri="http://doi.org/78.91011/j.celrep")],
                ["Consortia 1", "Consortia 2"],
            ),
            dataset_schema_version="3.1.0",
        )
        published_collection_3 = self.generate_published_collection(
            owner="other owner",
            curator_name="other curator",
            add_datasets=1,
            metadata=CollectionMetadata(
                "test_collection_3",
                "described",
                "john doe",
                "john.doe@email.com",
                [Link(name="doi link", type=CollectionLinkType.DOI.name, uri="http://doi.org/78.91011/j.celrep")],
                ["Consortia 1", "Consortia 2"],
            ),
            dataset_schema_version="4.0.0",
        )
        self.generate_unpublished_collection(add_datasets=4)
        self.generate_revision(published_collection_1.collection_id)

        sorted_dataset_ids = (
            [published_collection_3.datasets[0].dataset_id.id]
            + [published_collection_2.datasets[0].dataset_id.id]
            + sorted([d.dataset_id.id for d in published_collection_1.datasets], reverse=True)
        )

        with self.subTest("Query without parameter"):
            response = self.app.get("/curation/v1/datasets")
            received_dataset_ids = []
            for dataset in response.json:
                received_dataset_ids.append(dataset["dataset_id"])
            self.assertTrue(200, response.status_code)
            self.assertEqual(sorted_dataset_ids, received_dataset_ids)

        with self.subTest("Query against major schema version 3"):
            response = self.app.get("/curation/v1/datasets?schema_version=3")
            received_dataset_ids = []
            for dataset in response.json:
                received_dataset_ids.append(dataset["dataset_id"])
            self.assertTrue(200, response.status_code)
            self.assertEqual(sorted_dataset_ids[1:], received_dataset_ids)

        with self.subTest("Query against major schema version 4"):
            response = self.app.get("/curation/v1/datasets?schema_version=4")
            received_dataset_ids = []
            for dataset in response.json:
                received_dataset_ids.append(dataset["dataset_id"])
            self.assertTrue(200, response.status_code)
            self.assertEqual([published_collection_3.datasets[0].dataset_id.id], received_dataset_ids)

        with self.subTest("Query against minor schema version"):
            response = self.app.get("/curation/v1/datasets?schema_version=3.0")
            received_dataset_ids = []
            for dataset in response.json:
                received_dataset_ids.append(dataset["dataset_id"])
            self.assertTrue(200, response.status_code)
            self.assertEqual(sorted_dataset_ids[2:], received_dataset_ids)

        with self.subTest("Query against patch schema version"):
            response = self.app.get("/curation/v1/datasets?schema_version=3.1.0")
            received_dataset_ids = []
            for dataset in response.json:
                received_dataset_ids.append(dataset["dataset_id"])
            self.assertTrue(200, response.status_code)
            self.assertEqual([published_collection_2.datasets[0].dataset_id.id], received_dataset_ids)

        with self.subTest("Query against non-matching major schema version"):
            response = self.app.get("/curation/v1/datasets?schema_version=5")
            self.assertTrue(200, response.status_code)
            self.assertCountEqual(response.json, [])

        with self.subTest("Query against non-matching minor schema version"):
            response = self.app.get("/curation/v1/datasets?schema_version=3.2")
            self.assertTrue(200, response.status_code)
            self.assertCountEqual(response.json, [])

        with self.subTest("Query against non-matching patch schema version"):
            response = self.app.get("/curation/v1/datasets?schema_version=3.1.1")
            self.assertTrue(200, response.status_code)
            self.assertCountEqual(response.json, [])

        with self.subTest("Query against non-version input"):
            response = self.app.get("/curation/v1/datasets?schema_version=invalid_input")
            self.assertEqual(400, response.status_code)


class TestGetDatasetVersion(BaseAPIPortalTest):
    @patch("backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn)
    def test_get_dataset_version_ok(self, mock_config: Mock):
        collection = self.generate_published_collection()
        collection_id = collection.collection_id
        initial_published_dataset = collection.datasets[0]
        initial_published_dataset_version_id = collection.datasets[0].version_id
        published_revision = self.generate_revision(collection_id)
        published_dataset_revision = self.generate_dataset(
            collection_version=published_revision,
            replace_dataset_version_id=initial_published_dataset_version_id,
            publish=True,
        )

        headers = self.make_owner_header()
        # get previously published dataset version
        test_url = f"/curation/v1/dataset_versions/{initial_published_dataset_version_id}"
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(initial_published_dataset_version_id.id, response.json["dataset_version_id"])
        self.assertEqual(initial_published_dataset.dataset_id.id, response.json["dataset_id"])
        self.assertEqual(collection_id.id, response.json["collection_id"])
        expected_assets = [  # Filter out disallowed file types + properly construct url
            {
                "filesize": -1,
                "filetype": "H5AD",
                "url": f"http://domain/{initial_published_dataset_version_id}.h5ad",
            },
            {
                "filesize": -1,
                "filetype": "RDS",
                "url": f"http://domain/{initial_published_dataset_version_id}.rds",
            },
        ]
        self.assertEqual(response.json["assets"], expected_assets)

        # get currently published dataset version
        test_url = f"/curation/v1/dataset_versions/{published_dataset_revision.dataset_version_id}"
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(published_dataset_revision.dataset_version_id, response.json["dataset_version_id"])
        self.assertEqual(published_dataset_revision.dataset_id, response.json["dataset_id"])
        self.assertEqual(collection_id.id, response.json["collection_id"])
        expected_assets = [  # Filter out disallowed file types + properly construct url
            {
                "filesize": -1,
                "filetype": "H5AD",
                "url": f"http://domain/{published_dataset_revision.dataset_version_id}.h5ad",
            },
            {
                "filesize": -1,
                "filetype": "RDS",
                "url": f"http://domain/{published_dataset_revision.dataset_version_id}.rds",
            },
        ]
        self.assertEqual(response.json["assets"], expected_assets)

    def test_get_dataset_version_4xx(self):
        headers = self.make_owner_header()
        with self.subTest("Input is not valid UUID"):
            test_url = "/curation/v1/dataset_versions/bad-input-id"
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(403, response.status_code)
        with self.subTest("Input is valid (canonical dataset ID) but not found"):
            dataset = self.generate_published_collection().datasets[0]
            # passing in canonical dataset ID instead of version ID
            test_url = f"/curation/v1/dataset_versions/{dataset.dataset_id}"
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(404, response.status_code)
        with self.subTest("Input is ID for unpublished dataset version"):
            dataset = self.generate_unpublished_collection(add_datasets=1).datasets[0]
            test_url = f"/curation/v1/dataset_versions/{dataset.dataset_id}"
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(404, response.status_code)
        with self.subTest("Input is ID for unpublished dataset revision"):
            collection = self.generate_published_collection(add_datasets=1)
            unpublished_revision = self.generate_revision(collection.collection_id)
            unpublished_dataset_revision = self.generate_dataset(
                collection_version=unpublished_revision,
                replace_dataset_version_id=DatasetVersionId(collection.datasets[0].version_id),
            )
            test_url = f"/curation/v1/dataset_versions/{unpublished_dataset_revision.dataset_version_id}"
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(404, response.status_code)
        with self.subTest("Dataset is part of a tombstoned Collection"):
            collection = self.generate_published_collection()
            dataset = collection.datasets[0]
            self.business_logic.tombstone_collection(collection.collection_id)
            test_url = f"/curation/v1/dataset_versions/{dataset.version_id}"
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(410, response.status_code)
        with self.subTest("Dataset version is prior published; Collection is tombstoned"):
            collection = self.generate_published_collection(add_datasets=2)
            dataset_version_id = collection.datasets[0].version_id
            revision = self.generate_revision(collection.collection_id)
            self.business_logic.publish_collection_version(revision.version_id)
            self.business_logic.tombstone_collection(collection.collection_id)
            test_url = f"/curation/v1/dataset_versions/{dataset_version_id}"
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(410, response.status_code)


class TestGetDatasetIdVersions(BaseAPIPortalTest):
    def test_get_dataset_id_versions_ok(self):
        collection = self.generate_published_collection()
        collection_id = collection.collection_id
        dataset_id = collection.datasets[0].dataset_id
        dataset_version_id = collection.datasets[0].version_id
        published_revision = self.generate_revision(collection_id)
        published_dataset_revision = self.generate_dataset(
            collection_version=published_revision, replace_dataset_version_id=dataset_version_id, publish=True
        )
        unpublished_revision = self.generate_revision(collection_id)
        self.generate_dataset(
            collection_version=unpublished_revision,
            replace_dataset_version_id=DatasetVersionId(published_dataset_revision.dataset_version_id),
        )

        test_url = f"/curation/v1/datasets/{dataset_id}/versions"
        headers = self.make_owner_header()
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(200, response.status_code)
        expected = defaultdict(list)
        for dataset in response.json:
            self.assertIsNone(dataset.get("revised_at"))
            expected["dataset_version_ids"].append(dataset["dataset_version_id"])
            expected["collection_ids"].append(dataset["collection_id"])
            expected["collection_version_ids"].append(dataset["collection_version_id"])
            expected["published_at"].append(dataset["published_at"])
        # Check that only published datasets appear
        # Must be returned in reverse chronological order
        self.assertEqual(
            [published_dataset_revision.dataset_version_id, dataset_version_id.id], expected["dataset_version_ids"]
        )
        self.assertEqual([collection_id.id, collection_id.id], expected["collection_ids"])
        self.assertEqual(
            [published_revision.version_id.id, collection.version_id.id], expected["collection_version_ids"]
        )
        self.assertTrue(expected["published_at"][0] > expected["published_at"][1])

    def test_get_dataset_id_version_4xx(self):
        with self.subTest("Input is not a UUID"):
            test_url = "/curation/v1/datasets/not-uuid-input/versions"
            headers = self.make_owner_header()
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(403, response.status_code)
        with self.subTest("Dataset with that UUID does not exist"):
            test_url = f"/curation/v1/datasets/{str(uuid.uuid4())}/versions"
            headers = self.make_owner_header()
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(404, response.status_code)
        with self.subTest("Dataset Exists, but not published so no public version history"):
            dataset_id = self.generate_unpublished_collection(add_datasets=1).datasets[0].dataset_id
            test_url = f"/curation/v1/datasets/{dataset_id}/versions"
            headers = self.make_owner_header()
            response = self.app.get(test_url, headers=headers)
            self.assertEqual(404, response.status_code)
        with self.subTest("Dataset is part of a tombstoned Collection"):
            collection = self.generate_published_collection()
            dataset = collection.datasets[0]
            self.business_logic.tombstone_collection(collection.collection_id)
            test_url = f"/curation/v1/datasets/{dataset.dataset_id}/versions"
            response = self.app.get(test_url, headers=self.make_owner_header())
            self.assertEqual(410, response.status_code)
        with self.subTest("Dataset version is prior published; Collection is tombstoned"):
            collection = self.generate_published_collection(add_datasets=2)
            dataset_id = collection.datasets[0].dataset_id
            revision = self.generate_revision(collection.collection_id)
            self.business_logic.publish_collection_version(revision.version_id)
            self.business_logic.tombstone_collection(collection.collection_id)
            test_url = f"/curation/v1/datasets/{dataset_id}/versions"
            response = self.app.get(test_url, headers=self.make_owner_header())
            self.assertEqual(410, response.status_code)
        with self.subTest("Dataset is tombstoned"):
            collection = self.generate_published_collection(add_datasets=2)
            dataset = collection.datasets[0]
            revision = self.generate_revision(collection.collection_id)
            self.business_logic.remove_dataset_version(revision.version_id, dataset.version_id, delete_published=True)
            self.business_logic.publish_collection_version(revision.version_id)
            test_url = f"/curation/v1/datasets/{dataset.dataset_id}/versions"
            response = self.app.get(test_url, headers=self.make_owner_header())
            self.assertEqual(410, response.status_code)


class TestPostDataset(BaseAPIPortalTest):
    """
    Unit test for POST /datasets, which is used to add an empty dataset to a collection version
    """

    def test_post_datasets_nonexistent_collection_403(self):
        with self.subTest("UUID input valid, but not found"):
            non_existent_collection_id = str(uuid.uuid4())
            test_url = f"/curation/v1/collections/{non_existent_collection_id}/datasets"
            headers = self.make_owner_header()
            response = self.app.post(test_url, headers=headers)
            self.assertEqual(404, response.status_code)
        with self.subTest("UUID input invalid"):
            non_existent_collection_id = "123-example-fake-uuid"
            test_url = f"/curation/v1/collections/{non_existent_collection_id}/datasets"
            headers = self.make_owner_header()
            response = self.app.post(test_url, headers=headers)
            self.assertEqual(403, response.status_code)

    def test_post_datasets_with_collection_201(self):
        collection = self.generate_unpublished_collection()
        test_id = collection.collection_id
        test_url = f"/curation/v1/collections/{test_id}/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)
        self.assertTrue(response.json["dataset_id"])

    def test_post_datasets_with_collection_version_id_403(self):
        collection = self.generate_unpublished_collection()
        test_id = collection.version_id
        test_url = f"/curation/v1/collections/{test_id}/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)

    def test_post_datasets_super(self):
        collection = self.generate_unpublished_collection()
        test_id = collection.collection_id
        test_url = f"/curation/v1/collections/{test_id}/datasets"
        headers = self.make_super_curator_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)

    def test_post_datasets_not_owner_403(self):
        collection = self.generate_collection()
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets"
        headers = self.make_not_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)

    def test_post_datasets_public_collection_405(self):
        collection = self.generate_collection(visibility="PUBLIC")
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(405, response.status_code)

    def test_post_datasets_no_auth_401(self):
        collection = self.generate_collection(visibility="PUBLIC")
        test_url = f"/curation/v1/collections/{collection.version_id}/datasets"
        response = self.app.post(test_url)
        self.assertEqual(401, response.status_code)

    def test_post_datasets_returns_canonical_id(self):
        """
        POST /datasets returns the canonical dataset id on creation.
        """
        collection = self.generate_unpublished_collection()
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)

        looked_up_version = self.business_logic.get_collection_version(collection.version_id)
        self.assertEqual(1, len(looked_up_version.datasets))
        self.assertEqual(response.json["dataset_id"], looked_up_version.datasets[0].dataset_id.id)


class TestPostRevision(BaseAPIPortalTest):
    def test__post_revision__no_auth(self):
        collection_id = self.generate_collection(visibility=CollectionVisibility.PUBLIC.name).collection_id
        response = self.app.post(f"/curation/v1/collections/{collection_id}/revision")
        self.assertEqual(401, response.status_code)

    def test__post_revision__Not_Public(self):
        collection_id = self.generate_unpublished_collection().collection_id
        headers = self.make_super_curator_header()
        response = self.app.post(f"/curation/v1/collections/{collection_id.id}/revision", headers=headers)
        self.assertEqual(403, response.status_code)

    def test__post_revision__Not_Owner(self):
        collection_id = self.generate_collection(
            visibility=CollectionVisibility.PUBLIC.name, owner="someone else"
        ).collection_id
        response = self.app.post(
            f"/curation/v1/collections/{collection_id}/revision",
            headers=self.make_owner_header(),
        )
        self.assertEqual(403, response.status_code)

    def test__post_revision__OK(self):
        collection_id = self.generate_collection(visibility=CollectionVisibility.PUBLIC.name).collection_id
        response = self.app.post(
            f"/curation/v1/collections/{collection_id}/revision",
            headers=self.make_owner_header(),
        )
        self.assertEqual(201, response.status_code)
        self.assertNotEqual(collection_id, response.json["collection_id"])

    def test__post_revision_by_collection_version_id_403(self):
        version_id = self.generate_collection(visibility=CollectionVisibility.PUBLIC.name).version_id
        response = self.app.post(
            f"/curation/v1/collections/{version_id}/revision",
            headers=self.make_owner_header(),
        )
        self.assertEqual(403, response.status_code)

    def test__post_revision__Super_Curator(self):
        collection_id = self.generate_collection(visibility=CollectionVisibility.PUBLIC.name).collection_id
        headers = self.make_super_curator_header()
        response = self.app.post(f"/curation/v1/collections/{collection_id}/revision", headers=headers)
        self.assertEqual(201, response.status_code)
        self.assertNotEqual(collection_id, response.json["collection_id"])


@patch(
    "backend.common.utils.dl_sources.url.DropBoxURL.file_info",
    return_value={"size": 1, "name": "file.h5ad"},
)
@patch("backend.layers.thirdparty.step_function_provider.StepFunctionProvider")
class TestPutLink(BaseAPIPortalTest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        cls.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"

    def test__from_link__no_auth(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should fail with 401 Unauthorized if the user is not authenticated
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)]
        )
        body = {"link": self.good_link}
        headers = None
        for id in [dataset.dataset_version_id, dataset.dataset_id]:
            response = self.app.put(
                f"/curation/v1/collections/{dataset.collection_id}/datasets/{id}",
                json=body,
                headers=headers,
            )

            self.assertEqual(401, response.status_code)

    def test__from_link__Not_Public(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should fail with 403 Unauthorized
        if trying to upload to a published collection
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)],
            publish=True,
        )
        body = {"link": self.good_link}
        headers = self.make_owner_header()
        for id in [dataset.dataset_version_id, dataset.dataset_id]:
            response = self.app.put(
                f"/curation/v1/collections/{dataset.collection_id}/datasets/{id}",
                json=body,
                headers=headers,
            )

            self.assertEqual(403, response.status_code)

    def test__from_link__Not_Owner(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should fail with 403 Unauthorized
        if the authenticated user is not the owner of a collection
        """

        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)],
        )
        body = {"link": self.dummy_link}
        headers = self.make_not_owner_header()
        for id in [dataset.dataset_version_id, dataset.dataset_id]:
            response = self.app.put(
                f"/curation/v1/collections/{dataset.collection_id}/datasets/{id}",
                json=body,
                headers=headers,
            )

            self.assertEqual(403, response.status_code)

    def test__new_from_link__OK(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should succeed if a valid link is uploaded by the owner of the collection or
        a super curator
        """

        def _test_create(collection_id, dataset_id, headers):
            body = {"link": self.good_link}
            response = self.app.put(
                f"/curation/v1/collections/{collection_id}/datasets/{dataset_id}",
                json=body,
                headers=headers,
            )
            self.assertEqual(202, response.status_code)
            headers = [("owner", self.make_owner_header()), ("super curator", self.make_super_curator_header())]
            for auth_type, header in headers:
                with self.subTest(f"{auth_type}, unpublished collection"):
                    dataset = self.generate_dataset(
                        statuses=[
                            DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)
                        ],
                    )
                    _test_create(dataset.collection_id, dataset.dataset_id, header)

                with self.subTest(f"{auth_type}, revision"):
                    collection_id = self.generate_published_collection().collection_id
                    revision = self.generate_revision(collection_id)
                    dataset = self.generate_dataset(
                        statuses=[
                            DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)
                        ],
                        collection_version=revision,
                    )
                    _test_create(revision.version_id, dataset.dataset_id, header)

    def test__existing_from_link__OK(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id on an existing dataset_id should succeed if a valid link is uploaded by the
        owner of the collection or a super curator
        """

        def _test_create(collection_id, dataset_id, headers):
            body = {"link": self.good_link}
            response = self.app.put(
                f"/curation/v1/collections/{collection_id}/datasets/{dataset_id}",
                json=body,
                headers=headers,
            )
            self.assertEqual(202, response.status_code)
            headers = [("owner", self.make_owner_header()), ("super curator", self.make_super_curator_header())]
            for auth_type, header in headers:
                with self.subTest(f"{auth_type}, unpublished collection"):
                    dataset = self.generate_dataset(
                        statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
                    )
                    _test_create(dataset.collection_id, dataset.dataset_id, header)

                with self.subTest(f"{auth_type}, revision"):
                    collection_id = self.generate_published_collection().collection_id
                    revision = self.generate_revision(collection_id)
                    dataset = self.generate_dataset(
                        statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
                        collection_version=revision,
                    )
                    _test_create(revision.version_id, dataset.dataset_id, header)

    def test_from_link__403(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id  should fail if version IDs are used for unpublished collections or
        for referencing datasets
        """

        def _test_create(collection_version_id, dataset_version_id):
            body = {"link": self.good_link}
            headers = self.make_owner_header()
            response = self.app.put(
                f"/curation/v1/collections/{collection_version_id}/datasets/{dataset_version_id}",
                json=body,
                headers=headers,
            )
            self.assertEqual(403, response.status_code)

        with self.subTest("use collection version ID for unpublished collection"):
            collection = self.generate_unpublished_collection()
            dataset = self.generate_dataset(
                statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
                collection_version=collection,
            )
            _test_create(collection.version_id, dataset.dataset_id)

        with self.subTest("use dataset version ID for unpublished collection"):
            dataset = self.generate_dataset(
                statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
            )
            _test_create(dataset.collection_id, dataset.dataset_version_id)


class TestAuthToken(BaseAPIPortalTest):
    @patch("backend.curation.api.v1.curation.auth.token.CorporaAuthConfig")
    @patch("backend.curation.api.v1.curation.auth.token.auth0_management_session")
    def test__post_token__201(self, auth0_management_session: Mock, CorporaAuthConfig: Mock):
        test_secret = "password1234"
        test_email = "user@email.com"
        test_user_id = "test_user_id"
        CorporaAuthConfig().api_key_secret = test_secret
        auth0_management_session.get_user_api_key_identity = Mock(return_value={"profileData": {"email": test_email}})
        auth0_management_session.generate_access_token = Mock(return_value={"access_token": "OK"})
        user_api_key = generate(test_user_id, test_secret)
        response = self.app.post("/curation/v1/auth/token", headers={"x-api-key": user_api_key})
        self.assertEqual(201, response.status_code)
        token = response.json["access_token"]
        self.assertEqual("OK", token)
        auth0_management_session.get_user_api_key_identity.assert_called_once_with(test_user_id)

    @patch("backend.curation.api.v1.curation.auth.token.CorporaAuthConfig")
    def test__post_token__401(self, CorporaAuthConfig):
        test_secret = "password1234"
        test_user_id = "test_user_id"
        CorporaAuthConfig().api_key_secret = test_secret
        user_api_key = generate(test_user_id, "not the right secret")
        response = self.app.post("/curation/v1/auth/token", headers={"x-api-key": user_api_key})
        self.assertEqual(401, response.status_code)

    @patch("backend.curation.api.v1.curation.auth.token.CorporaAuthConfig")
    @patch("backend.curation.api.v1.curation.auth.token.auth0_management_session")
    def test__post_token__404(self, auth0_management_session, CorporaAuthConfig):
        test_secret = "password1234"
        test_user_id = "test_user_id"
        CorporaAuthConfig().api_key_secret = test_secret
        auth0_management_session.get_user_api_key_identity = Mock(return_value=None)
        user_api_key = generate(test_user_id, test_secret)
        response = self.app.post("/curation/v1/auth/token", headers={"x-api-key": user_api_key})
        self.assertEqual(404, response.status_code)
