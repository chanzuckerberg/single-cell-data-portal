import json
from unittest.mock import Mock, patch

from backend.common.corpora_orm import (
    CollectionVisibility,
    DbDataset,
    IsPrimaryData,
    ProcessingStatus,
    UploadStatus,
    ValidationStatus,
)
from backend.common.providers.crossref_provider import CrossrefDOINotFoundException
from backend.common.utils.api_key import generate
from backend.portal.api.curation.v1.curation.collections.common import EntityColumns
from unit.backend.api_server.base_api_test import BaseAuthAPITest
from unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from unit.backend.layers.common.base_api_test import NewBaseTest


class TestAsset(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()
        self.test_dataset_id = "test_dataset_id"

    def test__get_dataset_asset__OK(self):
        bucket = self.CORPORA_TEST_CONFIG["bucket_name"]
        s3_file_name = "test_s3_uri.h5ad"
        content = "Hello world!"
        self.create_s3_object(s3_file_name, bucket, content=content)

        expected_body = [dict(filename="test_filename", filesize=len(content), filetype="H5AD")]

        response = self.app.get(
            f"/curation/v1/collections/test_collection_id/datasets/{self.test_dataset_id}/assets",
        )
        self.assertEqual(200, response.status_code)
        actual_body = response.json
        presign_url = actual_body[0].pop("presigned_url")
        self.assertIsNotNone(presign_url)
        self.assertEqual(expected_body, actual_body)

    def test__get_dataset_asset__file_error(self):
        expected_body = [dict(filename="test_filename", filesize=-1, filetype="H5AD")]

        response = self.app.get(
            f"/curation/v1/collections/test_collection_id/datasets/{self.test_dataset_id}/assets",
        )
        self.assertEqual(202, response.status_code)
        actual_body = response.json
        presign_url = actual_body[0].pop("presigned_url", None)
        self.assertIsNone(presign_url)
        self.assertEqual(expected_body, actual_body)

    def test__get_dataset_asset__dataset_NOT_FOUND(self):
        bad_id = "bad_id"

        response = self.app.get(
            f"/curation/v1/collections/test_collection_id/datasets/{bad_id}/assets",
        )
        self.assertEqual(404, response.status_code)
        actual_body = response.json
        self.assertEqual("Dataset not found.", actual_body["detail"])

    def test__get_dataset_asset__collection_dataset_NOT_FOUND(self):
        """Return Not found when the dataset is not part of the collection"""
        bad_id = "bad_id"
        test_url = f"/curation/v1/collections/test_collection_id/datasets/{bad_id}/assets"
        response = self.app.get(test_url)
        self.assertEqual(404, response.status_code)
        actual_body = response.json
        self.assertEqual("Dataset not found.", actual_body["detail"])

    @patch(
        "backend.portal.api.curation.v1.curation.collections.collection_id.datasets.dataset_id."
        "assets.get_dataset_else_error",
        return_value=None,
    )
    def test__get_dataset_asset__asset_NOT_FOUND(self, get_dataset_else_error: Mock):
        mocked_dataset = Mock()
        mocked_dataset.get_assets.return_value = None
        get_dataset_else_error.return_value = mocked_dataset
        response = self.app.get(f"/curation/v1/collections/test_collection_id/datasets/{self.test_dataset_id}/assets")
        self.assertEqual(404, response.status_code)
        actual_body = response.json
        self.assertEqual("No assets found. The dataset may still be processing.", actual_body["detail"])
        mocked_dataset.get_assets.assert_called()


class TestDeleteCollection(BaseAuthAPITest):
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

    def test__delete_revision_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 204), ("super", 204)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                public_collection_id = self.generate_published_collection().collection_id
                revision_collection_id = self.generate_collection(
                    visibility=CollectionVisibility.PRIVATE.name, revision_of=public_collection_id
                ).collection_id
                self._test(revision_collection_id, auth, expected_response)

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
                tombstone_collection_id = self.generate_collection(
                    visibility=CollectionVisibility.PUBLIC.name, tombstone=True
                ).collection_id
                self._test(tombstone_collection_id, auth, expected_response)


class TestS3Credentials(NewBaseTest):
    @patch("backend.portal.api.curation.v1.curation.collections.collection_id.s3_upload_credentials.sts_client")
    def test__generate_s3_credentials__OK(self, sts_client: Mock):
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
            collection_id = self.generate_unpublished_collection().collection_id
            headers = {"Authorization": f"Bearer {token}"}

            response = self.app.get(f"/curation/v1/collections/{collection_id}/s3-upload-credentials", headers=headers)
            self.assertEqual(200, response.status_code)
            token_sub = self._mock_assert_authorized_token(token)["sub"]
            self.assertEqual(response.json["Bucket"], "cellxgene-dataset-submissions-test")
            if is_super_curator:
                self.assertEqual(response.json["UploadKeyPrefix"], f"super/{collection_id}/")
            else:
                self.assertEqual(response.json["UploadKeyPrefix"], f"{token_sub}/{collection_id}/")

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


class TestPostCollection(NewBaseTest):
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
        self.assertIn("id", response.json.keys())
        self.assertEqual(201, response.status_code)

    def test__create_collection__InvalidParameters(self):
        requests = [
            (
                dict(
                    name="",
                    description="",
                    contact_name="",
                    contact_email="@email.com",
                    doi="10.111/not_curie_reference_format",
                ),
                [
                    {"name": "contact_email", "reason": "Invalid format."},
                    {"name": "description", "reason": "Cannot be blank."},
                    {"name": "name", "reason": "Cannot be blank."},
                    {"name": "contact_name", "reason": "Cannot be blank."},
                    {"name": "DOI", "reason": "DOI must be a CURIE reference."},
                ],
                5,
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


class TestGetCollections(NewBaseTest):
    def setUp(self):
        super().setUp()
        self.test_collection = dict(
            name="collection", description="description", contact_name="john doe", contact_email="johndoe@email.com"
        )

    def test__get_collections_no_auth__OK(self):
        res_no_auth = self.app.get("/curation/v1/collections")
        self.assertEqual(200, res_no_auth.status_code)
        self.assertEqual(6, len(res_no_auth.json))
        [self.assertEqual("PUBLIC", c["visibility"]) for c in res_no_auth.json]

    def test__get_collections_with_auth__OK(self):
        res_auth = self.app.get("/curation/v1/collections", headers=self.make_owner_header())
        with self.subTest("The 'revising_in' attribute is None for unauthorized public collections"):
            conditions_tested = 0
            for c in res_auth.json:
                if c["id"] in (
                    "test_collection_id_public_for_revision_one",
                    "test_collection_id_public_for_revision_two",
                ):
                    conditions_tested += 1
                    self.assertIsNone(c["revising_in"])
            self.assertEqual(2, conditions_tested)
        with self.subTest("The 'revising_in' attribute is None for collections which lack a revision"):
            conditions_tested = 0
            for c in res_auth.json:
                if c["id"] in (
                    "test_collection_id_public",
                    "test_collection_with_link",
                    "test_collection_with_link_and_dataset_changes",
                ):
                    conditions_tested += 1
                    self.assertIsNone(c["revising_in"])
            self.assertEqual(3, conditions_tested)
        with self.subTest("The 'revising_in' attribute is equal to the id of the revision Collection"):
            conditions_tested = 0
            for c in res_auth.json:
                if c["id"] == "test_collection_id":
                    conditions_tested += 1
                    self.assertEqual("test_collection_id_revision", c["revising_in"])
            self.assertTrue(1, conditions_tested)
        with self.subTest("Datasets in a revision contain a reference to their published counterpart, if it exists"):
            conditions_tested = 0
            for c in res_auth.json:
                if c["id"] == "test_collection_id_revision":
                    for d in c["datasets"]:
                        if d["id"] == "test_publish_revision_with_links__revision_dataset":
                            conditions_tested += 1
                            self.assertEqual("test_dataset_id", d["revision_of"])
            self.assertTrue(1, conditions_tested)

    def test__get_collections_no_auth_visibility_private__OK(self):
        params = {"visibility": "PRIVATE"}
        res_private = self.app.get("/curation/v1/collections", query_string=params)
        self.assertEqual(403, res_private.status_code)

    def test__get_collections_no_auth_visibility_public__OK(self):
        params = {"visibility": "PUBLIC"}
        res_public = self.app.get("/curation/v1/collections", query_string=params)
        self.assertEqual(200, res_public.status_code)
        self.assertEqual(6, len(res_public.json))

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

    def test__get_only_public_collections_with_auth__OK(self):
        params = {"visibility": "PUBLIC"}
        res = self.app.get("/curation/v1/collections", query_string=params, headers=self.make_owner_header())
        self.assertEqual(200, res.status_code)
        self.assertEqual(6, len(res.json))
        [self.assertEqual("PUBLIC", c["visibility"]) for c in res.json]

    def test__get_only_private_collections_with_auth__OK(self):
        second_collection = self.generate_collection()
        for status in (ProcessingStatus.PENDING, ProcessingStatus.SUCCESS):
            self.generate_dataset(
                collection_id=second_collection.collection_id,
                processing_status={"processing_status": status},
            )
        third_collection = self.generate_collection()
        for status in (ProcessingStatus.INITIALIZED, ProcessingStatus.SUCCESS):
            self.generate_dataset(
                collection_id=third_collection.collection_id,
                processing_status={"processing_status": status},
            )
        params = {"visibility": "PRIVATE"}
        res = self.app.get("/curation/v1/collections", query_string=params, headers=self.make_owner_header())
        conditions_tested = 0
        with self.subTest("Summary collection-level processing statuses are accurate"):
            for collection in res.json:
                if collection["id"] in (second_collection.collection_id, third_collection.collection_id):
                    self.assertEqual(collection["processing_status"], "PENDING")
                    conditions_tested += 1
                else:
                    self.assertEqual(collection["processing_status"], "SUCCESS")
        self.assertEqual(2, conditions_tested)
        self.assertEqual(200, res.status_code)
        self.assertEqual(3, len(res.json))
        [self.assertEqual("PRIVATE", c["visibility"]) for c in res.json]

    def test__verify_expected_public_collection_fields(self):
        collection = self.generate_collection(
            visibility=CollectionVisibility.PUBLIC.name,
            links=[
                {
                    "link_name": "test_raw_data_link_name",
                    "link_type": "RAW_DATA",
                    "link_url": "http://test_raw_data_url.place",
                }
            ],
        )
        self.generate_dataset(collection_version_id=collection.version_id)
        res = self.app.get("/curation/v1/collections")
        self.assertEqual(200, res.status_code)
        for resp_collection in res.json:
            if resp_collection["id"] is collection.collection_id:
                break

        self.check_fields(EntityColumns.link_cols, resp_collection["links"][0], "links")
        self.check_fields(EntityColumns.dataset_preview_cols, resp_collection["datasets"][0], "datasets")

        collections_cols = EntityColumns.collections_cols.copy()
        collections_cols.remove("owner")
        collections_cols.remove("tombstone")
        collections_cols.append("processing_status")
        collections_cols.append("collection_url")
        collections_cols.append("doi")
        self.check_fields(collections_cols, resp_collection, "collection")

    def test__verify_expected_private_collection_fields(self):
        collection = self.generate_collection(
            visibility=CollectionVisibility.PRIVATE.name,
            links=[
                {
                    "link_name": "test_raw_data_link_name",
                    "link_type": "RAW_DATA",
                    "link_url": "http://test_raw_data_url.place",
                }
            ],
        )
        self.generate_dataset(collection_version_id=collection.version_id)
        params = {"visibility": "PRIVATE"}

        def _test(owner):
            if owner:
                header = self.make_owner_header()
                subtest_prefix = "owner"
            else:
                header = self.make_not_owner_header()
                subtest_prefix = "not_owner"
            res = self.app.get("/curation/v1/collections", query_string=params, headers=header)
            self.assertEqual(200, res.status_code)
            for resp_collection in res.json:
                if resp_collection["id"] is collection.collection_id:
                    break

            self.check_fields(EntityColumns.link_cols, resp_collection["links"][0], f"{subtest_prefix}:links")
            self.check_fields(
                EntityColumns.dataset_preview_cols, resp_collection["datasets"][0], f"{subtest_prefix}:datasets"
            )

            collections_cols = EntityColumns.collections_cols.copy()
            collections_cols.remove("tombstone")
            collections_cols.remove("owner")
            collections_cols.append("processing_status")
            collections_cols.append("collection_url")
            collections_cols.append("doi")
            self.check_fields(collections_cols, resp_collection, f"{subtest_prefix}:collection")

        _test(True)

    def check_fields(self, fields: list, response: dict, entity: str):
        for key in fields:
            with self.subTest(f"{entity}:{key}"):
                self.assertIn(key, response.keys())
                response.pop(key)
        with self.subTest(f"No Extra fields in {entity}"):
            self.assertFalse(response)

    def test__no_tombstoned_collections_or_datasets_included(self):
        second_collection = self.generate_collection(
            tombstone=False, name="second collection", visibility=CollectionVisibility.PUBLIC
        )
        self.generate_dataset(collection_id=second_collection.collection_id)
        self.generate_dataset(collection_id=second_collection.collection_id, tombstone=True)
        tombstoned_collection = self.generate_collection(
            tombstone=True, name="second collection", visibility=CollectionVisibility.PUBLIC
        )
        self.generate_dataset(collection_id=tombstoned_collection.collection_id, tombstone=True)

        res = self.app.get("/curation/v1/collections", headers=self.make_owner_header())

        contains_tombstoned_collection_flag = False
        for collection in res.json:
            if collection["id"] == second_collection.collection_id:
                self.assertEqual(1, len(collection["datasets"]))
            if collection["id"] == tombstoned_collection.collection_id:
                contains_tombstoned_collection_flag = True
        self.assertEqual(False, contains_tombstoned_collection_flag)


class TestGetCollectionID(NewBaseTest):
    expected_body = {
        "collection_url": "https://frontend.corporanet.local:3000/collections/test_collection_id",
        "contact_email": "somebody@chanzuckerberg.com",
        "contact_name": "Some Body",
        "curator_name": "",
        "datasets": [
            {
                "assay": [{"label": "test_assay", "ontology_term_id": "test_obo"}],
                "batch_condition": ["batchA", "batchB"],
                "cell_count": None,
                "cell_type": [{"label": "test_cell_type", "ontology_term_id": "test_opo"}],
                "dataset_assets": [{"filename": "test_filename", "filetype": "H5AD"}],
                "development_stage": [{"label": "test_development_stage", "ontology_term_id": "test_obo"}],
                "disease": [
                    {"label": "test_disease", "ontology_term_id": "test_obo"},
                    {"label": "test_disease2", "ontology_term_id": "test_obp"},
                    {"label": "test_disease3", "ontology_term_id": "test_obq"},
                ],
                "donor_id": ["donor_1", "donor_2"],
                "self_reported_ethnicity": [{"label": "test_ethnicity", "ontology_term_id": "test_obo"}],
                "explorer_url": "test_url",
                "id": "test_dataset_id",
                "is_primary_data": [True],
                "mean_genes_per_cell": 0.0,
                "title": "test_dataset_name",
                "organism": [{"label": "test_organism", "ontology_term_id": "test_obo"}],
                "processing_status": "PENDING",
                "revised_at": None,
                "revision": 0,
                "revision_of": None,
                "schema_version": "3.0.0",
                "sex": [
                    {"label": "test_sex", "ontology_term_id": "test_obo"},
                    {"label": "test_sex2", "ontology_term_id": "test_obp"},
                ],
                "suspension_type": ["nucleus"],
                "tissue": [{"label": "test_tissue", "ontology_term_id": "test_obo"}],
                "tombstone": False,
                "x_approximate_distribution": "NORMAL",
            }
        ],
        "description": "test_description",
        "doi": "",
        "id": "test_collection_id",
        "links": [
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
        ],
        "name": "test_collection_name",
        "processing_status": "PENDING",
        "published_at": None,
        "publisher_metadata": None,
        "revised_at": None,
        "revision_of": None,
        "revising_in": None,
        "tombstone": False,
        "visibility": "PUBLIC",
    }

    def setUp(self):
        super().setUp()
        self.test_collection = dict(
            name="collection", description="description", contact_name="john doe", contact_email="johndoe@email.com"
        )

    def test__get_public_collection_verify_body_is_reshaped_correctly__OK(self):
        dataset = self.session.query(DbDataset).filter(DbDataset.id == "test_dataset_id").one_or_none()
        self.assertIsInstance(dataset.organism, list)
        # Make this entry a dict instead of a list to test ability of the handler to reshape to list/array
        dataset.organism = dataset.organism[0]
        self.session.flush()
        dataset_modified = self.session.query(DbDataset).filter(DbDataset.id == "test_dataset_id").one_or_none()
        self.assertIsInstance(dataset_modified.organism, dict)

        res = self.app.get("/curation/v1/collections/test_collection_id")
        self.assertEqual(200, res.status_code)
        res_body = res.json
        del res_body["created_at"]  # too finicky; ignore
        self.assertDictEqual(self.expected_body, res_body)  # Confirm dict has been packaged in list
        self.assertEqual(json.dumps(self.expected_body, sort_keys=True), json.dumps(res_body))

    def test__get_private_collection__OK(self):
        res = self.app.get("/curation/v1/collections/test_collection_id_revision")
        self.assertEqual(200, res.status_code)
        self.assertEqual("test_collection_id_revision", res.json["id"])

    def test__get_collection_with_dataset_failing_validation(self):
        collection = self.generate_collection(
            visibility=CollectionVisibility.PRIVATE.name,
        )
        dataset = self.generate_dataset(
            collection=collection,
            processing_status={
                "processing_status": ProcessingStatus.FAILURE,
                "validation_status": ValidationStatus.INVALID,
                "validation_message": "test message",
            },
        )
        res = self.app.get(f"/curation/v1/collections/{collection.collection_id}")
        self.assertEqual("FAILURE", res.json["processing_status"])
        actual_dataset = res.json["datasets"][0]
        self.assertEqual(dataset.id, actual_dataset["id"])
        self.assertEqual("VALIDATION_FAILURE", actual_dataset["processing_status"])
        self.assertEqual("test message", actual_dataset["processing_status_detail"])

    def test__get_collection_with_dataset_failing_pipeline(self):
        collection = self.generate_collection(
            visibility=CollectionVisibility.PRIVATE.name,
        )
        dataset = self.generate_dataset(
            collection=collection, processing_status={"processing_status": ProcessingStatus.FAILURE}
        )
        res = self.app.get(f"/curation/v1/collections/{collection.collection_id}")
        self.assertEqual("FAILURE", res.json["processing_status"])
        actual_dataset = res.json["datasets"][0]
        self.assertEqual(dataset.id, actual_dataset["id"])
        self.assertEqual("PIPELINE_FAILURE", actual_dataset["processing_status"])

    def test__get_nonexistent_collection__403(self):
        res = self.app.get("/curation/v1/collections/test_collection_id_nonexistent")
        self.assertEqual(403, res.status_code)

    def test__get_tombstoned_collection__403(self):
        tombstoned_collection = self.generate_collection(
            tombstone=True, name="tombstoned collection", visibility=CollectionVisibility.PUBLIC
        )
        self.generate_dataset(collection_id=tombstoned_collection.collection_id, tombstone=True)
        res = self.app.get(f"/curation/v1/collections/{tombstoned_collection.collection_id}")
        self.assertEqual(403, res.status_code)

    def test_get_collection_with_no_datasets(self):
        collection = self.generate_collection(name="No Datasets", visibility=CollectionVisibility.PUBLIC)
        res = self.app.get(f"/curation/v1/collections/{collection.collection_id}")
        self.assertEqual(200, res.status_code)
        self.assertEqual(res.json["processing_status"], None)

    def test__get_collection_with_tombstoned_datasets__OK(self):
        collection = self.generate_collection(
            tombstone=False, name="collection", visibility=CollectionVisibility.PUBLIC
        )
        self.generate_dataset(collection_id=collection.collection_id, tombstone=False)
        self.generate_dataset(collection_id=collection.collection_id, tombstone=True)
        res = self.app.get(f"/curation/v1/collections/{collection.collection_id}")
        self.assertEqual(1, len(res.json["datasets"]))

    def test__get_public_collection_with_auth_access_type_write__OK(self):
        res = self.app.get("/curation/v1/collections/test_collection_id", headers=self.make_owner_header())
        self.assertEqual(200, res.status_code)
        self.assertEqual("test_collection_id", res.json["id"])

    def test__get_public_collection_with_auth_access_type_read__OK(self):
        res = self.app.get("/curation/v1/collections/test_collection_id_not_owner", headers=self.make_owner_header())
        self.assertEqual(200, res.status_code)
        self.assertEqual("test_collection_id_not_owner", res.json["id"])

    def test__get_private_collection_with_auth_access_type_write__OK(self):
        res = self.app.get("/curation/v1/collections/test_collection_id_revision", headers=self.make_owner_header())
        self.assertEqual(200, res.status_code)
        self.assertEqual("test_collection_id_revision", res.json["id"])

    def test__get_collectoin_with_x_approximate_distribution_none__OK(self):
        collection = self.generate_collection()
        self.generate_dataset(x_approximate_distribution=None, collection=collection)
        res = self.app.get(f"/curation/v1/collections/{collection.collection_id}", headers=self.make_owner_header())
        self.assertEqual(200, res.status_code)
        self.assertIsNone(res.json["datasets"][0]["x_approximate_distribution"])


class TestPatchCollectionID(NewBaseTest):
    def setUp(self):
        super().setUp()
        self.test_collection = dict(
            name="collection", description="description", contact_name="john doe", contact_email="johndoe@email.com"
        )

    def test__update_collection__no_auth(self):
        collection_id = self.generate_collection().collection_id
        response = self.app.patch(f"/curation/v1/collections/{collection_id}", data=json.dumps(self.test_collection))
        self.assertEqual(401, response.status_code)

    def test__update_collection__OK(self):
        collection_id = self.generate_collection().collection_id
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            data=json.dumps(self.test_collection),
            headers=self.make_owner_header(),
        )
        self.assertEqual(200, response.status_code)

    def test__update_collection_partial_data__OK(self):
        description = "a description"
        contact_name = "first last"
        contact_email = "name@server.domain"
        links = [
            {"link_name": "name", "link_type": "RAW_DATA", "link_url": "http://test_link.place"},
        ]
        publisher_metadata = {
            "authors": [{"name": "First Last", "given": "First", "family": "Last"}],
            "journal": "Journal of Anamolous Results",
            "published_year": 2022,
            "published_month": 1,
        }
        collection_id = self.generate_collection(
            description=description,
            contact_name=contact_name,
            contact_email=contact_email,
            links=links,
            publisher_metadata=publisher_metadata,
        ).collection_id
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
        self.assertEqual(description, response.json["description"])
        self.assertEqual(contact_name, response.json["contact_name"])
        self.assertEqual(contact_email, response.json["contact_email"])
        self.assertEqual(links, response.json["links"])
        self.assertEqual(publisher_metadata, response.json["publisher_metadata"])

    def test_update_collection_with_empty_required_fields(self):
        tests = [dict(description=""), dict(contact_name=""), dict(contact_email=""), dict(name="")]

        collection_id = self.generate_collection().collection_id
        for test in tests:
            with self.subTest(test):
                response = self.app.patch(
                    f"/curation/v1/collections/{collection_id}",
                    data=json.dumps(test),
                    headers=self.make_owner_header(),
                )
                self.assertEqual(400, response.status_code)

    def test__update_collection__links__OK(self):
        name = "partial updates test collection"
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
                collection_id = self.generate_collection(name=name, links=initial_links).collection_id
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
                    self.assertEqual(name, response.json["name"])
                    self.assertEqual(expected_links, response.json["links"])

    def test__update_collection__doi__OK(self):
        initial_doi = "12.3456/doi_curie_reference"
        links = [
            {"link_name": "new doi", "link_type": "DOI", "link_url": initial_doi},
        ]
        new_doi = "10.1016"  # a real DOI (CURIE reference)
        collection_id = self.generate_collection(links=links).collection_id
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

    def test__update_collection__doi_is_not_CURIE_reference__BAD_REQUEST(self):
        links = [
            {"link_name": "doi", "link_type": "DOI", "link_url": "http://doi.doi/10.1011/something"},
        ]
        collection = self.generate_collection(links=links)
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

    def test__update_collection__doi_does_not_exist__BAD_REQUEST(self):
        name = "bad doi update test collection"
        links = [
            {"link_name": "name", "link_type": "RAW_DATA", "link_url": "http://test_link.place"},
            {"link_name": "doi", "link_type": "DOI", "link_url": "http://doi.doi/10.1011/something"},
        ]
        new_links = [
            {"link_name": "new link", "link_type": "RAW_DATA", "link_url": "http://brand_new_link.place"},
        ]
        publisher_metadata = {
            "authors": [{"name": "First Last", "given": "First", "family": "Last"}],
            "journal": "Journal of Anamolous Results",
            "published_year": 2022,
            "published_month": 1,
        }
        with patch(
            "backend.common.providers.crossref_provider.CrossrefProvider.fetch_metadata",
            side_effect=CrossrefDOINotFoundException(),
        ):
            collection = self.generate_collection(name=name, links=links, publisher_metadata=publisher_metadata)
            collection_id = collection.collection_id
            original_collection = self.app.get(f"curation/v1/collections/{collection_id}").json

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
            self.assertEqual(name, original_collection_unchanged["name"])
            # Only compare to first item in links list because "DOI" type gets removed from Curator API response
            self.assertEqual(links[:1], original_collection_unchanged["links"])
            self.assertEqual(publisher_metadata, original_collection_unchanged["publisher_metadata"])

    def test__update_collection__Not_Owner(self):
        collection_id = self.generate_collection(owner="someone else").collection_id
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}",
            data=json.dumps(self.test_collection),
            headers=self.make_owner_header(),
        )
        self.assertEqual(403, response.status_code)

    def test__update_collection__Super_Curator(self):
        collection_id = self.generate_collection().collection_id
        headers = self.make_super_curator_header()
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}", data=json.dumps(self.test_collection), headers=headers
        )
        self.assertEqual(200, response.status_code)

    def test__update_public_collection_owner__405(self):
        collection_id = self.generate_collection(visibility=CollectionVisibility.PUBLIC).collection_id
        headers = self.make_super_curator_header()
        response = self.app.patch(
            f"/curation/v1/collections/{collection_id}", data=json.dumps(self.test_collection), headers=headers
        )
        self.assertEqual(405, response.status_code)
        self.assertEqual(
            "Directly editing a public Collection is not allowed; you must create a revision.",
            json.loads(response.text)["detail"],
        )


class TestDeleteDataset(BaseAuthAPITest):
    def test__delete_dataset(self):
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        auth_credentials = [
            (self.make_super_curator_header, "super", 202),
            (self.make_owner_header, "owner", 202),
            (None, "none", 401),
            (self.make_not_owner_header, "not_owner", 403),
        ]
        for auth, auth_description, expected_status_code in auth_credentials:
            with self.subTest(f"{auth_description} {expected_status_code}"):
                collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name)
                dataset = self.generate_dataset(
                    collection=collection,
                    processing_status=processing_status,
                )
                test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.id}"
                headers = auth() if callable(auth) else auth
                response = self.app.delete(test_url, headers=headers)
                self.assertEqual(expected_status_code, response.status_code)


class TestGetDatasets(BaseAuthAPITest):
    def test_get_dataset_in_a_collection_200(self):
        collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(collection_version_id=collection.version_id, name="test")
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.id}"

        response = self.app.get(test_url)
        self.assertEqual(200, response.status_code)
        self.assertEqual(dataset.id, response.json["id"])

    def test_get_dataset_shape(self):
        collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(collection_version_id=collection.version_id, name="test")
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.id}"
        response = self.app.get(test_url)
        self.assertEqual(dataset.name, response.json["title"])

    def test_get_dataset_is_primary_data_shape(self):
        collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name)
        tests = [
            (IsPrimaryData.PRIMARY, [True]),
            (IsPrimaryData.SECONDARY, [False]),
            (IsPrimaryData.BOTH, [True, False]),
        ]
        for is_primary_data, result in tests:
            with self.subTest(f"{is_primary_data}=={result}"):
                dataset = self.generate_dataset(
                    collection_version_id=collection.version_id, is_primary_data=is_primary_data
                )
                test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.id}"
                response = self.app.get(test_url)
                self.assertEqual(result, response.json["is_primary_data"])

    def test_get_nonexistent_dataset_404(self):
        collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name)
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/1234-1243-2134-234-1342"
        response = self.app.get(test_url)
        self.assertEqual(404, response.status_code)

    def test_get_tombstoned_dataset_in_a_collection_404(self):
        collection = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(collection_version_id=collection.version_id, tombstone=True)
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.id}"
        response = self.app.get(test_url)
        self.assertEqual(404, response.status_code)

    def test_get_datasets_nonexistent_collection_404(self):
        test_url = "/curation/v1/collections/nonexistent/datasets/1234-1243-2134-234-1342"
        headers = self.make_owner_header()
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(404, response.status_code)


class TestPostDataset(BaseAuthAPITest):
    def test_post_datasets_nonexistent_collection_403(self):
        test_url = "/curation/v1/collections/nonexistent/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)

    def test_post_datasets_201(self):
        collection = self.generate_collection()
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)
        self.assertTrue(response.json["id"])

    def test_post_datasets_super(self):
        collection = self.generate_collection()
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets"
        headers = self.make_super_curator_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)

    def test_post_datasets_not_owner_201(self):
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
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets"
        response = self.app.post(test_url)
        self.assertEqual(401, response.status_code)


class TestPostRevision(BaseAuthAPITest):
    def test__post_revision__no_auth(self):
        collection_id = self.generate_collection(visibility=CollectionVisibility.PUBLIC.name).collection_id
        response = self.app.post(f"/curation/v1/collections/{collection_id}/revision")
        self.assertEqual(401, response.status_code)

    def test__post_revision__Not_Public(self):
        collection_id = self.generate_collection().collection_id
        headers = self.make_super_curator_header()
        response = self.app.post(f"/curation/v1/collections/{collection_id}/revision", headers=headers)
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
        self.assertNotEqual(collection_id, response.json["id"])

    def test__post_revision__Super_Curator(self):
        collection_id = self.generate_collection(visibility=CollectionVisibility.PUBLIC.name).collection_id
        headers = self.make_super_curator_header()
        response = self.app.post(f"/curation/v1/collections/{collection_id}/revision", headers=headers)
        self.assertEqual(201, response.status_code)
        self.assertNotEqual(collection_id, response.json["id"])


@patch(
    "backend.common.utils.dl_sources.url.DropBoxURL.file_info",
    return_value={"size": 1, "name": "file.h5ad"},
)
@patch("backend.common.upload.start_upload_sfn")
class TestPutLink(BaseAuthAPITest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        cls.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"

    def _test_new(self, headers: dict = None, collection_params: dict = None, body: dict = None):
        headers = headers if headers else {}
        collection_params = collection_params if collection_params else {}
        collection = self.generate_collection(**collection_params)
        dataset = self.generate_dataset(
            collection_id=collection.collection_id,
            processing_status={"processing_status": ProcessingStatus.INITIALIZED},
        )
        body = body if body else ""
        headers["Content-Type"] = "application/json"
        response = self.app.put(
            f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.id}", json=body, headers=headers
        )
        return response

    def test__from_link__no_auth(self, *mocks):
        response = self._test_new(body={"link": self.good_link})
        self.assertEqual(401, response.status_code)

    def test__from_link__Not_Public(self, *mocks):
        headers = self.make_owner_header()
        response = self._test_new(
            headers, collection_params={"visibility": CollectionVisibility.PUBLIC}, body={"link": self.good_link}
        )
        self.assertEqual(403, response.status_code)

    def test__from_link__Not_Owner(self, *mocks):
        response = self._test_new(self.make_not_owner_header(), body={"link": self.dummy_link})
        self.assertEqual(403, response.status_code)

    def test__new_from_link__OK(self, *mocks):
        headers = self.make_owner_header()
        response = self._test_new(headers, body={"link": self.good_link})
        self.assertEqual(202, response.status_code)

    def test__new_from_link__Super_Curator(self, *mocks):
        headers = self.make_super_curator_header()
        response = self._test_new(headers, body={"link": self.good_link})
        self.assertEqual(202, response.status_code)

    def _test_existing(self, headers: dict = None):
        headers = headers if headers else {}
        headers["Content-Type"] = "application/json"
        collection = self.generate_collection()
        processing_status = dict(processing_status=ProcessingStatus.SUCCESS)
        dataset = self.generate_dataset(collection_id=collection.collection_id, processing_status=processing_status)
        body = {"link": self.good_link}
        response = self.app.put(
            f"/curation/v1/collections/{collection.collection_id}/datasets/{dataset.id}",
            data=json.dumps(body),
            headers=headers,
        )
        return response

    def test__existing_from_link__OK(self, *mocks):
        with self.subTest("dataset_id"):
            response = self._test_existing(self.make_owner_header())
            self.assertEqual(202, response.status_code)

    def test__existing_from_link__Super_Curator(self, *mocks):
        headers = self.make_super_curator_header()
        with self.subTest("dataset_id"):
            response = self._test_existing(headers)
            self.assertEqual(202, response.status_code)


class TestAuthToken(NewBaseTest):
    @patch("backend.portal.api.curation.v1.curation.auth.token.CorporaAuthConfig")
    @patch("backend.portal.api.curation.v1.curation.auth.token.auth0_management_session")
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

    @patch("backend.portal.api.curation.v1.curation.auth.token.CorporaAuthConfig")
    def test__post_token__401(self, CorporaAuthConfig):
        test_secret = "password1234"
        test_user_id = "test_user_id"
        CorporaAuthConfig().api_key_secret = test_secret
        user_api_key = generate(test_user_id, "not the right secret")
        response = self.app.post("/curation/v1/auth/token", headers={"x-api-key": user_api_key})
        self.assertEqual(401, response.status_code)

    @patch("backend.portal.api.curation.v1.curation.auth.token.CorporaAuthConfig")
    @patch("backend.portal.api.curation.v1.curation.auth.token.auth0_management_session")
    def test__post_token__404(self, auth0_management_session, CorporaAuthConfig):
        test_secret = "password1234"
        test_user_id = "test_user_id"
        CorporaAuthConfig().api_key_secret = test_secret
        auth0_management_session.get_user_api_key_identity = Mock(return_value=None)
        user_api_key = generate(test_user_id, test_secret)
        response = self.app.post("/curation/v1/auth/token", headers={"x-api-key": user_api_key})
        self.assertEqual(404, response.status_code)
