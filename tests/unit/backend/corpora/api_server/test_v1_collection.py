import itertools
import json
from datetime import datetime
from furl import furl
import unittest

from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
    UploadStatus,
    generate_uuid,
)
from backend.corpora.common.entities import Collection
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token


class TestCollection(BaseAuthAPITest):
    def validate_collections_response_structure(self, body):
        self.assertIn("collections", body)
        self.assertTrue(all(k in ["collections", "from_date", "to_date"] for k in body))

        for collection in body["collections"]:
            self.assertListEqual(sorted(collection.keys()), ["created_at", "id", "visibility"])
            self.assertEqual(collection["visibility"], "PUBLIC")
            self.assertGreaterEqual(datetime.fromtimestamp(collection["created_at"]).year, 1969)

    def validate_collection_uuid_response_structure(self, body):
        required_keys = [
            "name",
            "description",
            "id",
            "visibility",
            "links",
            "datasets",
            "created_at",
            "updated_at",
            "obfuscated_uuid",
            "contact_email",
            "contact_name",
            "data_submission_policy_version",
            "access_type",
        ]
        self.assertListEqual(sorted(body.keys()), sorted(required_keys))
        self.assertGreaterEqual(datetime.fromtimestamp(body["created_at"]).year, 1969)
        self.assertGreaterEqual(datetime.fromtimestamp(body["updated_at"]).year, 1969)

        for link in body["links"]:
            self.assertListEqual(sorted(link.keys()), ["link_name", "link_type", "link_url"])

        for dataset in body["datasets"]:
            required_keys = [
                "id",
                "assay",
                "tissue",
                "disease",
                "sex",
                "ethnicity",
                "organism",
                "development_stage",
                "name",
                "revision",
                "dataset_deployments",
                "dataset_assets",
                "processing_status",
                "created_at",
                "updated_at",
                "collection_id",
                "collection_visibility",
                "is_valid",
                "cell_count",
            ]
            self.assertListEqual(sorted(dataset.keys()), sorted(required_keys))

    def test__list_collection_options__allow(self):
        origin = "http://localhost:3000"
        res = self.app.options("/dp/v1/collections", headers={"origin": origin})
        self.assertEqual(res.status_code, 200)
        self.assertEqual(origin, res.headers["Access-Control-Allow-Origin"])

    def test__list_collection_options__no_allow(self):
        res = self.app.options("/dp/v1/collections", headers={"origin": "http://localhost:ABCD"})
        self.assertEqual(res.status_code, 200)
        self.assertNotIn("Access-Control-Allow-Origin", res.headers.keys())

    def test__list_collection__ok(self):
        path = "/dp/v1/collections"
        headers = dict(host="localhost")

        from_date = int(datetime.fromtimestamp(60).timestamp())
        creation_time = 70
        to_date = int(datetime.fromtimestamp(80).timestamp())

        expected_id = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, created_at=datetime.fromtimestamp(creation_time)
        ).id

        with self.subTest("No Parameters"):
            test_url = furl(path=path)
            response = self.app.get(test_url.url, headers=headers)
            self.assertEqual(200, response.status_code)
            actual_body = json.loads(response.data)
            self.validate_collections_response_structure(actual_body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["collections"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(None, actual_body.get("from_date"))

        with self.subTest("from_date"):
            test_url = furl(path=path, query_params={"from_date": from_date})
            response = self.app.get(test_url.url, headers=headers)
            self.assertEqual(200, response.status_code)
            actual_body = json.loads(response.data)
            self.validate_collections_response_structure(actual_body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["collections"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(actual_body["from_date"], from_date)

        with self.subTest("to_date"):
            test_url = furl(path=path, query_params={"to_date": to_date})
            response = self.app.get(test_url.url, headers=headers)
            self.assertEqual(200, response.status_code)
            actual_body = json.loads(response.data)
            self.validate_collections_response_structure(actual_body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["collections"]])
            self.assertEqual(to_date, actual_body["to_date"])
            self.assertEqual(None, actual_body.get("from_date"))

        with self.subTest("from_date->to_date"):
            test_url = furl(path=path, query_params={"from_date": from_date, "to_date": to_date})
            response = self.app.get(test_url.url, headers=headers)
            self.assertEqual(200, response.status_code)
            actual_body = json.loads(response.data)
            self.validate_collections_response_structure(actual_body)
            self.assertEqual(expected_id, actual_body["collections"][0]["id"])
            self.assertEqual(creation_time, actual_body["collections"][0]["created_at"])
            self.assertEqual(from_date, actual_body["from_date"])
            self.assertEqual(to_date, actual_body["to_date"])

    def test__get_collection_uuid__ok(self):
        """Verify the test collection exists and the expected fields exist."""
        expected_body = {
            "datasets": [
                {
                    "assay": [{"ontology_term_id": "test_obo", "label": "test_assay"}],
                    "dataset_assets": [
                        {
                            "dataset_id": "test_dataset_id",
                            "filename": "test_filename",
                            "filetype": "H5AD",
                            "id": "test_dataset_artifact_id",
                            "s3_uri": "s3://bogus-bucket/test_s3_uri.h5ad",
                            "type": "ORIGINAL",
                            "user_submitted": True,
                        }
                    ],
                    "dataset_deployments": [{"url": "test_url"}],
                    "development_stage": [{"label": "test_development_stage", "ontology_term_id": "test_obo"}],
                    "disease": [
                        {"label": "test_disease", "ontology_term_id": "test_obo"},
                        {"label": "test_disease2", "ontology_term_id": "test_obp"},
                        {"label": "test_disease3", "ontology_term_id": "test_obq"},
                    ],
                    "ethnicity": [{"label": "test_ethnicity", "ontology_term_id": "test_obo"}],
                    "linked_genesets": ["test_geneset_with_dataset"],
                    "id": "test_dataset_id",
                    "is_primary_data": "PRIMARY",
                    "mean_genes_per_cell": 0.0,
                    "name": "test_dataset_name",
                    "organism": {"label": "test_organism", "ontology_term_id": "test_obo"},
                    "collection_id": "test_collection_id",
                    "collection_visibility": "PUBLIC",
                    "cell_type": [{"label": "test_cell_type", "ontology_term_id": "test_opo"}],
                    "x_normalization": "test_x_normalization",
                    "x_approximate_distribution": "NORMAL",
                    "is_valid": False,
                    "revision": 0,
                    "sex": [
                        {"label": "test_sex", "ontology_term_id": "test_obo"},
                        {"label": "test_sex2", "ontology_term_id": "test_obp"},
                    ],
                    "tissue": [{"label": "test_tissue", "ontology_term_id": "test_obo"}],
                    "tombstone": False,
                    "processing_status": {
                        "id": "test_dataset_processing_status_id",
                        "dataset_id": "test_dataset_id",
                        "processing_status": "PENDING",
                        "upload_status": "UPLOADING",
                        "upload_progress": 4 / 9,
                        "validation_status": "NA",
                        "conversion_loom_status": "NA",
                        "conversion_anndata_status": "NA",
                        "conversion_rds_status": "NA",
                        "conversion_cxg_status": "NA",
                    },
                    "published": False,
                    "tombstone": False,
                    "schema_version": "2.0.0",
                }
            ],
            "description": "test_description",
            "genesets": [
                {
                    "collection_id": "test_collection_id",
                    "collection_visibility": "PUBLIC",
                    "linked_datasets": [],
                    "description": "this is a geneset",
                    "id": "test_geneset",
                    "name": "test_geneset",
                },
                {
                    "collection_id": "test_collection_id",
                    "collection_visibility": "PUBLIC",
                    "linked_datasets": ["test_dataset_id"],
                    "description": "this is a geneset with a dataset",
                    "id": "test_geneset_with_dataset",
                    "name": "test_geneset_with_dataset",
                },
            ],
            "id": "test_collection_id",
            "links": [
                {"link_name": "test_doi_link_name", "link_type": "DOI", "link_url": "http://test_doi_url.place"},
                {"link_name": "", "link_type": "DOI", "link_url": "http://test_no_link_name_doi_url.place"},
                {
                    "link_name": "test_raw_data_link_name",
                    "link_type": "RAW_DATA",
                    "link_url": "http://test_raw_data_url.place",
                },
                {"link_name": "", "link_type": "RAW_DATA", "link_url": "http://test_no_link_name_raw_data_url.place"},
                {
                    "link_name": "test_protocol_link_name",
                    "link_type": "PROTOCOL",
                    "link_url": "http://test_protocol_url.place",
                },
                {"link_name": "", "link_type": "PROTOCOL", "link_url": "http://test_no_link_name_protocol_url.place"},
                {
                    "link_name": "test_lab_website_link_name",
                    "link_type": "LAB_WEBSITE",
                    "link_url": "http://test_lab_website_url.place",
                },
                {
                    "link_name": "",
                    "link_type": "LAB_WEBSITE",
                    "link_url": "http://test_no_link_name_lab_website_url.place",
                },
                {"link_name": "test_other_link_name", "link_type": "OTHER", "link_url": "http://test_other_url.place"},
                {"link_name": "", "link_type": "OTHER", "link_url": "http://test_no_link_name_other_url.place"},
                {
                    "link_name": "test_data_source_link_name",
                    "link_type": "DATA_SOURCE",
                    "link_url": "http://test_data_source_url.place",
                },
                {
                    "link_name": "",
                    "link_type": "DATA_SOURCE",
                    "link_url": "http://test_no_link_name_data_source_url.place",
                },
            ],
            "name": "test_collection_name",
            "visibility": "PUBLIC",
            "obfuscated_uuid": "",
            "contact_name": "Some Body",
            "contact_email": "somebody@chanzuckerberg.com",
            "data_submission_policy_version": "0",
        }

        with self.subTest("auth cookie"):
            expected_body["access_type"] = "WRITE"
            test_url = furl(path="/dp/v1/collections/test_collection_id", query_params=dict(visibility="PUBLIC"))
            cxguser_cookie = get_auth_token(self.app)
            response = self.app.get(test_url.url, headers=dict(host="localhost", Cookie=cxguser_cookie))
            self.assertEqual(200, response.status_code)
            actual_body = self.remove_timestamps(json.loads(response.data))
            self.assertDictEqual(actual_body, expected_body)

        with self.subTest("no auth cookie"):
            expected_body["access_type"] = "READ"
            test_url = furl(path="/dp/v1/collections/test_collection_id", query_params=dict(visibility="PUBLIC"))
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            self.assertEqual(200, response.status_code)
            actual_body = self.remove_timestamps(json.loads(response.data))
            self.assertDictEqual(actual_body, expected_body)

    def test_get_collection_minimal__ok(self):
        with self.subTest("No Datasets"):
            collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name)
            test_url = furl(path=f"/dp/v1/collections/{collection.id}")
            resp = self.app.get(test_url.url)
            actual_body = self.remove_timestamps(json.loads(resp.data))
            expected_body = self.remove_timestamps(dict(**collection.reshape_for_api(), access_type="READ"))
            self.assertEqual(expected_body.pop("visibility").name, actual_body.pop("visibility"))
            self.assertEqual(expected_body, actual_body)

        with self.subTest("With a minimal dataset"):
            collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name)
            dataset = self.generate_dataset(
                self.session,
                collection_id=collection.id,
                organism=None,
                tissue=None,
                assay=None,
                disease=None,
                sex=None,
                ethnicity=None,
                development_stage=None,
                explorer_url=None,
            )
            expected_body = {
                "access_type": "READ",
                "contact_email": "",
                "contact_name": "",
                "data_submission_policy_version": "0",
                "datasets": [
                    {
                        "collection_id": collection.id,
                        "collection_visibility": "PUBLIC",
                        "dataset_assets": [],
                        "dataset_deployments": [],
                        "linked_genesets": [],
                        "name": dataset.name,
                        "id": dataset.id,
                        "is_valid": False,
                        "published": False,
                        "revision": 0,
                        "tombstone": False,
                    }
                ],
                "description": "",
                "genesets": [],
                "id": collection.id,
                "links": [],
                "name": "",
                "obfuscated_uuid": "",
                "visibility": "PUBLIC",
            }
            test_url = furl(path=f"/dp/v1/collections/{collection.id}")

            resp = self.app.get(test_url.url)

            self.assertEqual(200, resp.status_code)
            actual_body = self.remove_timestamps(json.loads(resp.data))
            expected_body = self.remove_timestamps(dict(**collection.reshape_for_api(), access_type="READ"))

            # Remove `visibility`, `collection_visibility` and `x_approximate_distribution` from the bodies,
            # since one will be an Enum and the other will be a string. That's expected since internal objects
            # should keep stricter data types, but JSON doesn't support Enums and therefore have to be converted
            # as strings.
            self.assertEqual(expected_body.pop("visibility").name, actual_body.pop("visibility"))
            self.assertEqual(
                expected_body["datasets"][0].pop("collection_visibility").name,
                actual_body["datasets"][0].pop("collection_visibility"),
            )
            self.assertEqual(
                expected_body["datasets"][0].pop("x_approximate_distribution").name,
                actual_body["datasets"][0].pop("x_approximate_distribution"),
            )
            self.assertEqual(
                expected_body["datasets"][0].pop("is_primary_data").name,
                actual_body["datasets"][0].pop("is_primary_data"),
            )

            self.assertEqual(expected_body, actual_body)

    def test__get_collection__ok(self):
        # Generate test cases
        authenticated = [True, False]
        owns = [True, False]
        visibility = ["public", "private"]
        obfuscated = [False]  # TODO: Once obfuscated uuid are support add True.
        test_cases = [params for params in itertools.product(authenticated, owns, visibility, obfuscated)]

        # Generate test collection
        test_collections = dict(
            public_not_owner=self.generate_collection(
                self.session, visibility=CollectionVisibility.PUBLIC.name, owner="someone else"
            ).id,
            public_owned=self.generate_collection(
                self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
            ).id,
            private_not_owner=self.generate_collection(
                self.session, visibility=CollectionVisibility.PRIVATE.name, owner="someone else"
            ).id,
            private_owned=self.generate_collection(
                self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
            ).id,
        )

        # run
        for auth, owns, visi, obfu in test_cases:
            expected_response_code = 200

            test_collection_id = test_collections["_".join([visi, "owned" if owns else "not_owner"])]
            expected_access_type = "WRITE" if owns and auth else "READ"

            with self.subTest(f"auth:{auth}, owns:{owns}, visi:{visi}, obfu:{obfu}, acc:{expected_access_type}"):
                if obfu:
                    raise NotImplementedError()
                test_url = furl(path=f"/dp/v1/collections/{test_collection_id}")
                if visi:
                    test_url.add(query_params=dict(visibility=visi.upper()))

                headers = dict(host="localhost")
                if auth:
                    headers["Cookie"] = get_auth_token(self.app)
                response = self.app.get(test_url.url, headers=headers)
                self.assertEqual(expected_response_code, response.status_code)
                if expected_response_code == 200:
                    actual_body = json.loads(response.data)
                    self.assertEqual(expected_access_type, actual_body["access_type"])

    def test_collection_with_tombstoned_dataset(self):
        dataset = self.generate_dataset(self.session, tombstone=True)
        test_url = furl(path="/dp/v1/collections/test_collection_id", query_params=dict(visibility="PUBLIC"))
        response = self.app.get(test_url.url, headers=dict(host="localhost"))

        self.assertEqual(response.status_code, 200)
        actual_dataset_ids = [d_id["id"] for d_id in json.loads(response.data)["datasets"]]

        self.assertNotIn(dataset.id, actual_dataset_ids)

    def test__get_collection_uuid__403_not_found(self):
        """Verify the test collection exists and the expected fields exist."""
        test_url = furl(path="/dp/v1/collections/AAAA-BBBB-CCCC-DDDD", query_params=dict(visibility="PUBLIC"))
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, response.status_code)

    def test__post_collection_returns_uuid_on_success(self):
        test_url = furl(path="/dp/v1/collections/")
        data = json.dumps(
            {
                "name": "collection name",
                "description": "This is a test collection",
                "contact_name": "person human",
                "contact_email": "person@human.com",
                "data_submission_policy_version": "0.0.1",
                "links": [
                    {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                    {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
                ],
            }
        )
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)},
            data=data,
        )
        self.assertEqual(201, response.status_code)

    def test__post_collection_fails_with_extra_fields(self):
        test_url = furl(path="/dp/v1/collections/")
        test_data = [
            {
                "name": "extra field in collection",
                "owner": "someone else",
                "description": "This is a test collection",
                "contact_name": "person human",
                "contact_email": "person@human.com",
                "data_submission_policy_version": "0.0.1",
                "links": [
                    {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                    {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
                ],
            },
            {
                "name": "extra field in link",
                "description": "This is a test collection",
                "contact_name": "person human",
                "contact_email": "person@human.com",
                "data_submission_policy_version": "0.0.1",
                "links": [
                    {
                        "link_name": "DOI Link",
                        "link_url": "http://doi.org/10.1017",
                        "link_type": "DOI",
                        "collection_id": "1235",
                    },
                ],
            },
        ]
        for data in test_data:
            with self.subTest(data["name"]):
                response = self.app.post(
                    test_url.url,
                    headers={
                        "host": "localhost",
                        "Content-Type": "application/json",
                        "Cookie": get_auth_token(self.app),
                    },
                    data=json.dumps(data),
                )
                self.assertEqual(400, response.status_code)

    def test__post_collection_fails_if_data_missing(self):
        test_url = furl(path="/dp/v1/collections/")
        data = json.dumps({"name": "bkjbjbjmbjm"})
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)},
            data=data,
        )
        self.assertEqual(400, response.status_code)

    def test__can_retrieve_created_collection(self):
        test_url = furl(path="/dp/v1/collections/")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        data = {
            "name": "another collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "data_submission_policy_version": "0.0.1",
            "links": [
                {"link_url": "http://doi.org/10.1016", "link_type": "OTHER"},
                {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
            ],
        }
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(data))
        self.assertEqual(201, response.status_code)
        collection_uuid = json.loads(response.data)["collection_uuid"]

        test_url = furl(path=f"/dp/v1/collections/{collection_uuid}")
        test_url.add(query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        self.assertEqual(body["description"], data["description"])
        self.assertEqual(body["name"], data["name"])
        self.assertEqual(body["contact_name"], body["contact_name"])
        self.assertEqual(body["contact_email"], body["contact_email"])

        # test that non owners only have read access
        no_cookie_headers = {"host": "localhost", "Content-Type": "application/json"}
        test_url = furl(path=f"/dp/v1/collections/{collection_uuid}")
        test_url.add(query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers=no_cookie_headers)
        self.assertEqual("READ", json.loads(response.data)["access_type"])

    def test__list_collection__check_owner(self):

        # Generate test collection
        public_owned = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        ).id
        private_owned = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        ).id
        public_not_owned = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="someone else"
        ).id
        private_not_owned = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="someone else"
        ).id

        path = "/dp/v1/collections"
        with self.subTest("no auth"):
            headers = {"host": "localhost", "Content-Type": "application/json"}
            response = self.app.get(path, headers=headers)
            self.assertEqual(200, response.status_code)
            result = json.loads(response.data)
            collections = result.get("collections")
            self.assertIsNotNone(collections)
            ids = [collection.get("id") for collection in collections]
            self.assertIn(public_owned, ids)
            self.assertIn(public_not_owned, ids)
            self.assertNotIn(private_owned, ids)
            self.assertNotIn(private_not_owned, ids)

        with self.subTest("auth"):
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.get(path, headers=headers)
            self.assertEqual(200, response.status_code)
            result = json.loads(response.data)
            collections = result.get("collections")
            self.assertIsNotNone(collections)
            ids = [collection.get("id") for collection in collections]
            self.assertIn(public_owned, ids)
            self.assertIn(public_not_owned, ids)
            self.assertIn(private_owned, ids)
            self.assertNotIn(private_not_owned, ids)


class TestCollectionDeletion(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def test_delete_private_collection__ok(self):
        # Generate test collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        processing_status_1 = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        processing_status_2 = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 100.0}

        dataset_1 = self.generate_dataset_with_s3_resources(
            self.session, collection=collection, processing_status=processing_status_1
        )
        dataset_2 = self.generate_dataset_with_s3_resources(
            self.session, collection=collection, processing_status=processing_status_2
        )

        s3_objects = self.get_s3_object_paths_from_dataset(dataset_1) + self.get_s3_object_paths_from_dataset(dataset_2)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        test_url = furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers=headers)

        self.assertEqual(response.status_code, 200)

        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset_1.id, dataset_ids)
        self.assertIn(dataset_2.id, dataset_ids)

        # delete collection
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 204)

        # check collection and datasets delete
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

        # check s3_resources have been deleted
        for bucket, key in s3_objects:
            self.assertS3FileDoesNotExist(bucket, key)

    @unittest.skip("Tombstone not supported through API.")
    def test_tombstone_collection__ok(self):
        # Generate test collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        # Generate the public collection with the same id as the private so a tombstone is created
        self.generate_collection(
            self.session, id=collection.id, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )

        processing_status_1 = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        processing_status_2 = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 100.0}

        dataset_1 = self.generate_dataset(self.session, collection=collection, processing_status=processing_status_1)
        dataset_2 = self.generate_dataset(self.session, collection=collection, processing_status=processing_status_2)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        test_url = furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PRIVATE"))

        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)

        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset_1.id, dataset_ids)
        self.assertIn(dataset_2.id, dataset_ids)

        # delete collection
        response = self.app.delete(test_url.url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # check collection and datasets delete
        self.session.expire_all()
        collection = Collection.get_collection(
            self.session, collection.id, CollectionVisibility.PRIVATE.name, include_tombstones=True
        )

        self.assertTrue(collection.tombstone)
        self.assertTrue(dataset_1.tombstone)
        self.assertTrue(dataset_2.tombstone)

        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_delete_collection__dataset_not_available(self):
        # Generate test collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        # Generate the public collection with the same id as the private so a tombstone is created
        self.generate_collection(
            self.session, id=collection.id, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 100.0}

        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        dataset_url = furl(path=f"/dp/v1/datasets/{dataset.id}/status")
        response = self.app.get(dataset_url.url, headers=headers)
        self.assertEqual(response.status_code, 200)

        test_url = furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.delete(test_url.url, headers=headers)

        self.assertEqual(response.status_code, 204)

        response = self.app.get(dataset_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_delete_collection__already_tombstoned__403(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id", tombstone=True
        )
        # Generate the public collection with the same id as the private so a tombstone is created
        self.generate_collection(
            self.session, id=collection.id, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )
        test_url = furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PRIVATE"))
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 204)

    def test_delete_collection__public__405(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )

        test_urls = [
            furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PUBLIC")),
            furl(path=f"/dp/v1/collections/{collection.id}"),
        ]
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        for test_url in test_urls:
            with self.subTest(test_url.url):
                response = self.app.delete(test_url.url, headers=headers)
                self.assertEqual(response.status_code, 405)

    def test_delete_collection__not_owner(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="someone_else"
        )
        test_url = furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PRIVATE"))
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_delete_collection__does_not_exist(self):
        fake_uuid = generate_uuid()
        test_url = furl(path=f"/dp/v1/collections/{fake_uuid}", query_params=dict(visibility="PRIVATE"))
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_deleted_collection_does_not_appear_in_collection_lists(self):
        private_collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        public_collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )
        collection_to_delete = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get("/dp/v1/collections/", headers=headers)

        collection_ids = [collection["id"] for collection in json.loads(response.data)["collections"]]
        self.assertIn(private_collection.id, collection_ids)
        self.assertIn(public_collection.id, collection_ids)
        self.assertIn(collection_to_delete.id, collection_ids)

        test_url = furl(path=f"/dp/v1/collections/{collection_to_delete.id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 204)

        # check not returned privately
        response = self.app.get("/dp/v1/collections/", headers=headers)
        collection_ids = [collection["id"] for collection in json.loads(response.data)["collections"]]
        self.assertIn(private_collection.id, collection_ids)
        self.assertIn(public_collection.id, collection_ids)

        self.assertNotIn(collection_to_delete.id, collection_ids)

        # check not returned publicly
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get("/dp/v1/collections/", headers=headers)
        collection_ids = [collection["id"] for collection in json.loads(response.data)["collections"]]
        self.assertIn(public_collection.id, collection_ids)
        self.assertNotIn(private_collection.id, collection_ids)
        self.assertNotIn(collection_to_delete.id, collection_ids)


class TestUpdateCollection(BaseAuthAPITest):
    def test__update_collection__OK(self):
        collection = self.generate_collection(self.session)
        test_fields = [
            "name",
            "description",
            "contact_name",
            "contact_email",
            "links",
            "data_submission_policy_version",
        ]
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        # Update the collection
        expected_body = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
            ],
            "data_submission_policy_version": "0",
        }
        data = json.dumps(expected_body)
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        for field in test_fields:
            self.assertEqual(expected_body[field], actual_body[field])

    def test__update_collection__403(self):
        collection = self.generate_collection(self.session, owner="someone else")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        data = json.dumps({"name": "new name"})
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)
        self.assertEqual(403, response.status_code)

    def test__update_collection_links__OK(self):
        collection = self.generate_collection(self.session)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        # add links
        links = [
            {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
            {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
        ]
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])

        # remove links
        links.pop()
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])

        # update links
        links = [{"link_name": "New name", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}]
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])

        # all together
        links = [
            {"link_name": "Link 1", "link_url": "This is a new link", "link_type": "OTHER"},
            {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
        ]
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])

        # Clear All links
        links = []
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)

        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])
