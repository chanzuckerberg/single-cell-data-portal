import copy
import json
import uuid
from dataclasses import asdict
from unittest.mock import Mock, patch

from backend.common.utils.api_key import generate
from backend.curation.api.v1.curation.collections.common import EntityColumns
from backend.layers.common.entities import (
    CollectionId,
    CollectionVersion,
    CollectionVisibility,
    DatasetArtifactType,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    Link,
    OntologyTermId,
)
from backend.layers.thirdparty.crossref_provider import CrossrefDOINotFoundException
from tests.unit.backend.layers.api.test_portal_api import generate_mock_publisher_metadata
from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest
from tests.unit.backend.layers.common.base_test import DatasetArtifactUpdate, DatasetStatusUpdate


class TestAsset(BaseAPIPortalTest):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()

    def test__get_dataset_asset__OK(self):
        self.business_logic.s3_provider.get_file_size = Mock(return_value=1000)
        self.business_logic.s3_provider.generate_presigned_url = Mock(return_value="http://mock.uri/presigned_url")

        dataset = self.generate_dataset(
            artifacts=[
                DatasetArtifactUpdate(DatasetArtifactType.H5AD, "http://mock.uri/asset.h5ad"),
                DatasetArtifactUpdate(DatasetArtifactType.CXG, "http://mock.uri/asset.cxg"),
                DatasetArtifactUpdate(DatasetArtifactType.RAW_H5AD, "http://mock.uri/raw.h5ad"),
            ],
            publish=True,
        )

        expected_body = [dict(filename="asset.h5ad", filesize=1000, filetype="H5AD")]

        response = self.app.get(
            f"/curation/v1/collections/test_collection_id/datasets/{dataset.dataset_id}/assets",
        )
        self.assertEqual(200, response.status_code)
        actual_body = response.json
        presign_url = actual_body[0].pop("presigned_url")
        self.assertIsNotNone(presign_url)
        self.assertEqual(expected_body, actual_body)

    def test__get_dataset_asset__file_error(self):
        self.business_logic.s3_provider.get_file_size = Mock(return_value=None)
        self.business_logic.s3_provider.generate_presigned_url = Mock(return_value=None)

        dataset = self.generate_dataset(
            artifacts=[
                DatasetArtifactUpdate(DatasetArtifactType.H5AD, "http://mock.uri/asset.h5ad"),
                DatasetArtifactUpdate(DatasetArtifactType.CXG, "http://mock.uri/asset.cxg"),
                DatasetArtifactUpdate(DatasetArtifactType.RAW_H5AD, "http://mock.uri/raw.h5ad"),
            ],
            publish=True,
        )

        expected_body = [dict(filename="asset.h5ad", filesize=-1, filetype="H5AD")]

        response = self.app.get(
            f"/curation/v1/collections/test_collection_id/datasets/{dataset.dataset_id}/assets",
        )
        self.assertEqual(202, response.status_code)
        actual_body = response.json
        presign_url = actual_body[0].pop("presigned_url", None)
        self.assertEqual(presign_url, "Not Found.")
        self.assertEqual(expected_body, actual_body)

    def test__get_dataset_asset__dataset_NOT_FOUND(self):
        collection = self.generate_published_collection()
        bad_id = uuid.uuid4()
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{bad_id}/assets"
        response = self.app.get(test_url)
        self.assertEqual(404, response.status_code)
        actual_body = response.json
        self.assertEqual("Dataset not found.", actual_body["detail"])

    def test__get_dataset_asset__asset_NOT_FOUND(self):
        dataset = self.generate_dataset(
            artifacts=[],
        )
        response = self.app.get(
            f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_id}/assets"
        )
        self.assertEqual(404, response.status_code)
        actual_body = response.json
        self.assertEqual("No assets found. The dataset may still be processing.", actual_body["detail"])


class TestDeleteCollection(BaseAPIPortalTest):
    def _test(self, collection_id, header, expected_status):
        if header == "owner":
            headers = self.make_owner_header()
        elif header == "super":
            headers = self.make_super_curator_header()
        elif header == "not_owner":
            headers = self.make_not_owner_header()
        elif "noauth":
            headers = {}

        response = self.app.delete(f"/curation/v1/collections/{collection_id}", headers=headers)
        self.assertEqual(expected_status, response.status_code)
        if response.status_code == 204:
            response = self.app.delete(f"/curation/v1/collections/{collection_id}", headers=headers)
            self.assertEqual(403, response.status_code)

    def test__delete_public_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 405), ("super", 405)]
        public_collection_id = self.generate_published_collection().collection_id
        for auth, expected_response in tests:
            with self.subTest(auth):
                self._test(public_collection_id, auth, expected_response)

    def test__delete_revision_collection_by_collection_version_id(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 204), ("super", 204)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                revision_collection = self.generate_collection_revision()
                self._test(revision_collection.version_id.id, auth, expected_response)

    def test__delete_private_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 204), ("super", 204)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                private_collection_id = self.generate_unpublished_collection().collection_id
                self._test(private_collection_id, auth, expected_response)

    def test__delete_tombstone_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 403), ("super", 403)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                collection_id = self.generate_published_collection().collection_id
                self.business_logic.tombstone_collection(collection_id)
                self._test(collection_id, auth, expected_response)


class TestS3Credentials(BaseAPIPortalTest):
    @patch("backend.common.corpora_config.CorporaConfig.__getattr__")
    @patch("backend.curation.api.v1.curation.collections.collection_id.s3_upload_credentials.sts_client")
    def test__generate_s3_credentials__OK(self, sts_client: Mock, mock_config: Mock):
        def mock_config_fn(name):
            if name == "curator_role_arn":
                return "test_role_arn"
            if name == "submission_bucket":
                return "cellxgene-dataset-submissions-test"

        mock_config.side_effect = mock_config_fn

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
            unpublished_collection = self.generate_unpublished_collection()
            headers = {"Authorization": f"Bearer {token}"}

            for id in [unpublished_collection.collection_id, unpublished_collection.version_id]:
                response = self.app.get(f"/curation/v1/collections/{id}/s3-upload-credentials", headers=headers)
                self.assertEqual(200, response.status_code)
                token_sub = self._mock_assert_authorized_token(token)["sub"]
                self.assertEqual(response.json["Bucket"], "cellxgene-dataset-submissions-test")
                if is_super_curator:
                    self.assertEqual(response.json["UploadKeyPrefix"], f"super/{id}/")
                else:
                    self.assertEqual(response.json["UploadKeyPrefix"], f"{token_sub}/{id}/")

        with self.subTest("collection owner"):
            _test("owner")

        with self.subTest("super curator"):
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


class TestGetCollectionID(BaseAPIPortalTest):
    def test__get_collection_verify_body_is_reshaped_correctly__OK(self):

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
            artifacts=[DatasetArtifactUpdate(type="h5ad", uri="http://test_filename")],
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
                "processing_status": "INITIALIZED",
                "tombstone": False,
                "processing_status_detail": None,
                "dataset_assets": [{"filename": "test_filename", "filetype": "H5AD"}],
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
                "processing_status": "PENDING",
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
        del res_body["datasets"][0]["revised_at"]  # too finicky; ignore

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

    def test_get_collection_dynamic_fields(self):
        def _test_responses_are_equal(test_ids, expected_id, auth_headers) -> dict:
            last_resp = None
            for header in auth_headers:
                """Check that a request with different permissions levels returns same response."""
                for identifier in test_ids:
                    """Check that a request with the listed IDs return the same response."""
                    resp = self.app.get(f"/curation/v1/collections/{identifier}", headers=header)
                    self.assertEqual(200, resp.status_code)
                    self.assertTrue(resp.json["collection_url"].endswith(expected_id.id))
                    if not last_resp:
                        last_resp = resp
                    else:
                        self.assertEqual(resp.json, last_resp.json, "All of the response bodies should be the same.")
            return last_resp.json

        privilage_access_headers = [self.make_owner_header(), self.make_super_curator_header()]
        restricted_access_headers = [self.make_not_owner_header(), self.make_not_auth_header()]
        all_headers = [*privilage_access_headers, *restricted_access_headers]
        unpublished = self.generate_unpublished_collection(add_datasets=1)
        with self.subTest("get unpublished version"):
            resp_collection = _test_responses_are_equal(
                [unpublished.version_id, unpublished.collection_id], unpublished.collection_id, all_headers
            )
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

        published = self.generate_published_collection(add_datasets=1)
        with self.subTest("get published version"):
            resp_collection = _test_responses_are_equal(
                [published.version_id, published.collection_id], published.collection_id, all_headers
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

        revision = self.generate_revision(published.collection_id)
        with self.subTest("get published with unpublished version and restricted access"):
            resp_collection = _test_responses_are_equal(
                [published.version_id, published.collection_id], published.collection_id, restricted_access_headers
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

        with self.subTest("get published with unpublished version and privilaged access"):
            resp_collection = _test_responses_are_equal(
                [published.version_id, published.collection_id], published.collection_id, privilage_access_headers
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

        with self.subTest("get unpublished version with published version and read access"):
            resp_collection = _test_responses_are_equal([revision.version_id], revision.version_id, all_headers)
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

        revised_dataset = self.generate_dataset(
            collection_version=revision, replace_dataset_version_id=revision.datasets[0].version_id
        )
        with self.subTest("get unpublished version with replace dataset and published version"):
            resp_collection = _test_responses_are_equal([revision.version_id], revision.version_id, all_headers)
            self.assertEqual("PRIVATE", resp_collection["visibility"])
            self.assertIsNone(resp_collection.get("revising_in"))
            self.assertEqual(revision.collection_id.id, resp_collection["revision_of"])
            self.assertIsNone(resp_collection["revised_at"])
            self.assertIsNotNone(resp_collection["published_at"])
            self.assertTrue(resp_collection["collection_url"].endswith(revision.version_id.id))
            self.assertEqual(revision.version_id.id, resp_collection["collection_version_id"])
            self.assertIsNotNone(resp_collection["datasets"][0]["revised_at"])
            self.assertEqual(revised_dataset.dataset_version_id, resp_collection["datasets"][0]["dataset_version_id"])
            self.assertIn(revised_dataset.dataset_version_id, resp_collection["datasets"][0]["explorer_url"])

        self.business_logic.publish_collection_version(revision.version_id)
        with self.subTest("get updated published version"):
            resp_collection = _test_responses_are_equal(
                [published.collection_id, revision.version_id], published.collection_id, all_headers
            )
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

    def test__get_nonexistent_collection__403(self):
        with self.subTest("UUID input valid, but not found"):
            non_existent_id = str(uuid.uuid4())
            res = self.app.get(f"/curation/v1/collections/{non_existent_id}")
            self.assertEqual(403, res.status_code)
        with self.subTest("UUID input invalid"):
            non_existent_id = "123-example-fake-uuid"
            res = self.app.get(f"/curation/v1/collections/{non_existent_id}")
            self.assertEqual(403, res.status_code)

    def test__get_tombstoned_collection__403(self):
        collection_version = self.generate_published_collection()
        self.business_logic.tombstone_collection(collection_version.collection_id)
        self._test_response(collection_version, 403)

    def test_get_collection_with_no_datasets(self):
        collection_version = self.generate_unpublished_collection(add_datasets=0)
        self._test_response(collection_version)

    def test_get_colletion_with_dataset_no_metadata(self):
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
        headers = self.make_owner_header() if auth else self.make_not_owner_header()
        res_1 = self.app.get(f"/curation/v1/collections/{version_id}", headers=headers)
        self.assertEqual(status_code, res_1.status_code)
        if status_code == 200:
            self.assertEqual(collection_version.collection_id.id, res_1.json["collection_id"])

        res_2 = self.app.get(f"/curation/v1/collections/{collection_id}", headers=headers)
        self.assertEqual(status_code, res_2.status_code)
        if status_code == 200:
            self.assertEqual(collection_version.collection_id.id, res_2.json["collection_id"])

        self.assertEqual(res_1.json, res_2.json)
        return res_1.json

    def test__get_collection_with_x_approximate_distribution_none__OK(self):
        metadata = copy.deepcopy(self.sample_dataset_metadata)
        metadata.x_approximate_distribution = None
        dataset = self.generate_dataset(metadata=metadata, publish=True)
        res = self.app.get(f"/curation/v1/collections/{dataset.collection_id}", headers=self.make_owner_header())
        self.assertEqual(200, res.status_code)
        self.assertIsNone(res.json["datasets"][0]["x_approximate_distribution"])


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

    def _delete(self, auth, collection_id, dataset_id):
        """
        Helper method to call the delete endpoint
        """
        test_url = f"/curation/v1/collections/{collection_id}/datasets/{dataset_id}"
        headers = auth() if callable(auth) else auth
        return self.app.delete(test_url, headers=headers)

    def test__delete_dataset_by_version_id(self):
        """
        Calling DELETE /collections/:collection_id/datasets/:dataset_id should work according to the
        auth token passed and when using versioned ids
        """
        for auth, auth_description, expected_status_code in self.auth_credentials:
            with self.subTest(f"{auth_description} {expected_status_code}"):
                dataset = self.generate_dataset(
                    statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)],
                    publish=False,
                )
                response = self._delete(auth, dataset.collection_version_id, dataset.dataset_version_id)
                self.assertEqual(expected_status_code, response.status_code)

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


class TestGetDatasets(BaseAPIPortalTest):
    def test_get_dataset_in_a_collection_200(self):
        dataset = self.generate_dataset(name="test")

        with self.subTest("by canonical collection_id"):
            test_url = f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}"

            response = self.app.get(test_url)
            self.assertEqual(200, response.status_code)
            self.assertEqual(dataset.dataset_id, response.json["dataset_id"])

        with self.subTest("by version_id"):
            test_url = f"/curation/v1/collections/{dataset.collection_version_id}/datasets/{dataset.dataset_version_id}"

            response = self.app.get(test_url)
            self.assertEqual(200, response.status_code)
            self.assertEqual(dataset.dataset_id, response.json["dataset_id"])

    def test_get_dataset_shape(self):
        # retrieve a private dataset
        private_dataset = self.generate_dataset(name="test")
        test_url = (
            f"/curation/v1/collections/{private_dataset.collection_id}/datasets/{private_dataset.dataset_version_id}"
        )
        response = self.app.get(test_url)
        body = response.json
        self.assertEqual("test", body["title"])

        # retrieve a public dataset
        public_dataset = self.generate_dataset(name="test", publish=True)
        test_url = (
            f"/curation/v1/collections/{public_dataset.collection_id}/datasets/{public_dataset.dataset_version_id}"
        )
        response = self.app.get(test_url)
        body = response.json
        self.assertEqual("test", body["title"])

        # retrieve a revised dataset using version_id
        collection_id = self.generate_published_collection(add_datasets=2).canonical_collection.id
        version = self.generate_revision(collection_id)
        dataset_version = self.generate_dataset(
            collection_version=version, replace_dataset_version_id=version.datasets[0].version_id
        )
        test_url = f"/curation/v1/collections/{version.version_id}/datasets/{dataset_version.dataset_version_id}"
        response = self.app.get(test_url)
        body = response.json

        # retrieve an unrevised dataset in a revision Collection
        unreplaced_dataset = version.datasets[1]
        test_url = f"/curation/v1/collections/{version.version_id}/datasets/{unreplaced_dataset.version_id}"
        response = self.app.get(test_url)
        body = response.json

        # retrieve a newly added dataset in a revision Collection
        new_dataset = self.generate_dataset(collection_version=version)
        test_url = f"/curation/v1/collections/{version.version_id}/datasets/{new_dataset.dataset_version_id}"
        response = self.app.get(test_url)
        body = response.json

        # retrieve a revision using dataset_id
        test_url = f"/curation/v1/collections/{version.version_id}/datasets/{dataset_version.dataset_id}"
        response = self.app.get(test_url)
        body = response.json

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
                test_url = f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}"
                response = self.app.get(test_url)
                self.assertEqual(result, response.json["is_primary_data"])

    def test_get_nonexistent_dataset_404(self):
        collection = self.generate_unpublished_collection()
        with self.subTest("UUID input valid, but not found"):
            non_existent_dataset_id = str(uuid.uuid4())
            test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{non_existent_dataset_id}"
            response = self.app.get(test_url)
            self.assertEqual(404, response.status_code)
        with self.subTest("UUID input invalid"):
            non_existent_dataset_id = "123-example-fake-uuid"
            test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{non_existent_dataset_id}"
            response = self.app.get(test_url)
            self.assertEqual(404, response.status_code)

    def test_get_datasets_nonexistent_collection_404(self):
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
            self.assertEqual(404, response.status_code)


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
            self.assertEqual(403, response.status_code)
        with self.subTest("UUID input invalid"):
            non_existent_collection_id = "123-example-fake-uuid"
            test_url = f"/curation/v1/collections/{non_existent_collection_id}/datasets"
            headers = self.make_owner_header()
            response = self.app.post(test_url, headers=headers)
            self.assertEqual(403, response.status_code)

    def test_post_datasets_with_collection_201(self):
        collection = self.generate_unpublished_collection()
        test_ids = [(collection.version_id, "version_id"), (collection.collection_id, "canonical_collection_id")]
        for test_id, test_name in test_ids:
            test_url = f"/curation/v1/collections/{test_id}/datasets"
            with self.subTest(test_name):
                headers = self.make_owner_header()
                response = self.app.post(test_url, headers=headers)
                self.assertEqual(201, response.status_code)
                self.assertTrue(response.json["dataset_id"])

    def test_post_datasets_super(self):
        collection = self.generate_unpublished_collection()
        test_ids = [(collection.version_id, "version_id"), (collection.collection_id, "canonical_collection_id")]
        for test_id, test_name in test_ids:
            test_url = f"/curation/v1/collections/{test_id}/datasets"
            with self.subTest(test_name):
                headers = self.make_super_curator_header()
                response = self.app.post(test_url, headers=headers)
                self.assertEqual(201, response.status_code)

    def test_post_datasets_not_owner_403(self):
        collection = self.generate_collection()
        test_url = f"/curation/v1/collections/{collection.version_id}/datasets"
        headers = self.make_not_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)

    def test_post_datasets_public_collection_405(self):
        collection = self.generate_collection(visibility="PUBLIC")
        test_url = f"/curation/v1/collections/{collection.version_id}/datasets"
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
        test_url = f"/curation/v1/collections/{collection.version_id}/datasets"
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
@patch("backend.common.upload.start_upload_sfn")
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
                f"/curation/v1/collections/{dataset.collection_version_id}/datasets/{id}",
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
                f"/curation/v1/collections/{dataset.collection_version_id}/datasets/{id}",
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
                f"/curation/v1/collections/{dataset.collection_version_id}/datasets/{id}",
                json=body,
                headers=headers,
            )

            self.assertEqual(403, response.status_code)

    def test__new_from_link__OK(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should succeed if a valid link is uploaded by the owner of the collection
        """

        def _test_create(collection_id, dataset_id):
            body = {"link": self.good_link}
            headers = self.make_owner_header()
            response = self.app.put(
                f"/curation/v1/collections/{collection_id}/datasets/{dataset_id}",
                json=body,
                headers=headers,
            )
            self.assertEqual(202, response.status_code)

        with self.subTest("with version_ids"):
            dataset = self.generate_dataset(
                statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)],
            )
            _test_create(dataset.collection_version_id, dataset.dataset_version_id)

        with self.subTest("with collection_ids"):
            dataset = self.generate_dataset(
                statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)],
            )
            _test_create(dataset.collection_id, dataset.dataset_id)

    def test__new_from_link__Super_Curator(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id should succeed if a valid link is uploaded by a super curator
        """

        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.INITIALIZED)],
        )
        body = {"link": self.good_link}
        headers = self.make_super_curator_header()
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_version_id}/datasets/{dataset.dataset_version_id}",
            json=body,
            headers=headers,
        )
        self.assertEqual(202, response.status_code)

    def test__existing_from_link__OK(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id on an existing dataset_id
        should succeed if a valid link is uploaded by the owner of the collection
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
        )
        body = {"link": self.good_link}
        headers = self.make_owner_header()
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_version_id}/datasets/{dataset.dataset_version_id}",
            json=body,
            headers=headers,
        )
        self.assertEqual(202, response.status_code)

    def test__existing_from_link__Super_Curator(self, *mocks):
        """
        Calling PUT /datasets/:dataset_id on an existing dataset_id
        should succeed if a valid link is uploaded by a super curator
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
        )
        body = {"link": self.good_link}
        headers = self.make_super_curator_header()
        response = self.app.put(
            f"/curation/v1/collections/{dataset.collection_version_id}/datasets/{dataset.dataset_version_id}",
            json=body,
            headers=headers,
        )
        self.assertEqual(202, response.status_code)


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
