import dataclasses
import itertools
import json
from typing import List
import unittest
from datetime import datetime
from unittest import mock
from unittest.mock import Mock, patch
from backend.layers.business.entities import DatasetArtifactDownloadData
from backend.layers.common.entities import CollectionId, CollectionVersionId, DatasetMetadata, DatasetProcessingStatus, DatasetUploadStatus, DatasetVersionId, Link, OntologyTermId
from backend.layers.thirdparty.uri_provider import FileInfo

from furl import furl

from backend.common.providers.crossref_provider import CrossrefDOINotFoundException, CrossrefFetchException
from backend.common.utils.corpora_constants import CorporaConstants
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from tests.unit.backend.layers.common.base_api_test import BaseAuthAPITest, DatasetArtifactUpdate, DatasetStatusUpdate, NewBaseTest


def generate_mock_publisher_metadata(journal_override=None):
    return {
        "authors": [{"given": "John", "family": "Doe"}, {"given": "Jane", "family": "Doe"}],
        "published_year": 2021,
        "published_month": 11,
        "published_day": 10,
        "published_at": 1636520400.0,
        "journal": journal_override or "Nature",
        "is_preprint": False,
    }


class TestCollection(NewBaseTest):

    # TODO: does not belong here
    def validate_collections_response_structure(self, body):
        self.assertIn("collections", body)
        self.assertTrue(all(k in ["collections", "from_date", "to_date"] for k in body))

        for collection in body["collections"]:
            self.assertListEqual(sorted(collection.keys()), ["created_at", "id", "visibility"])
            self.assertEqual(collection["visibility"], "PUBLIC")
            self.assertGreaterEqual(datetime.fromtimestamp(collection["created_at"]).year, 1969)

    # TODO: does not belong here
    def validate_collection_id_response_structure(self, body):
        required_keys = [
            "name",
            "description",
            "id",
            "visibility",
            "links",
            "datasets",
            "created_at",
            "updated_at",
            "contact_email",
            "contact_name",
            "curator_name",
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
                "self_reported_ethnicity",
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
                "is_valid",
                "cell_count",
            ]
            self.assertListEqual(sorted(dataset.keys()), sorted(required_keys))

    # ✅
    def test__list_collection_options__allow(self):
        origin = "http://localhost:3000"
        res = self.app.options("/dp/v1/collections", headers={"origin": origin})
        self.assertEqual(res.status_code, 200)
        self.assertEqual(origin, res.headers["Access-Control-Allow-Origin"])

    # ✅
    def test__list_collection_options__no_allow(self):
        res = self.app.options("/dp/v1/collections", headers={"origin": "http://localhost:ABCD"})
        self.assertEqual(res.status_code, 200)
        self.assertNotIn("Access-Control-Allow-Origin", res.headers.keys())

    # ✅
    def test__list_collection__ok(self):
        # TODO: anything to add now that the date filtering is deprecated?
        path = "/dp/v1/collections"
        headers = dict(host="localhost")

        # TODO: goes somewhere else

        expected_id = self.generate_published_collection().collection_id.id

        with self.subTest("No Parameters"):
            test_url = furl(path=path)
            response = self.app.get(test_url.url, headers=headers)
            self.assertEqual(200, response.status_code)
            actual_body = json.loads(response.data)
            # self.validate_collections_response_structure(actual_body) # TODO: maybe later
            self.assertIn(expected_id, [p["id"] for p in actual_body["collections"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(None, actual_body.get("from_date"))

    # ✅
    @patch("backend.layers.persistence.persistence_mock.datetime")
    def test__get_collection_id__ok(self, mock_dt):
        """Verify the test collection exists and the expected fields exist."""

        mock_dt.utcnow.return_value = 1234

        collection = self.generate_published_collection(add_datasets=2)

        self.maxDiff = None

        expected_body = {
            "access_type": "WRITE",
            "contact_email": "john.doe@email.com",
            "contact_name": "john doe",
            "created_at": 1234,
            "curator_name": "",
            "data_submission_policy_version": "1.0",
            "datasets": [
                {
                    "assay": [
                        {
                            "label": "test_assay_label",
                            "ontology_term_id": "test_assay_term_id"
                        }
                    ],
                    "batch_condition": [
                        "test_batch_1",
                        "test_batch_2"
                    ],
                    "cell_count": 10,
                    "cell_type": [
                        {
                            "label": "test_cell_type_label",
                            "ontology_term_id": "test_cell_type_term_id"
                        }
                    ],
                    "collection_id": collection.collection_id.id,
                    "created_at": 1234,
                    "dataset_assets": [],
                    "dataset_deployments": [
                        {
                            "url": "TODO"
                        }
                    ],
                    "development_stage": [
                        {
                            "label": "test_development_stage_label",
                            "ontology_term_id": "test_development_stage_term_id"
                        }
                    ],
                    "disease": [
                        {
                            "label": "test_disease_label",
                            "ontology_term_id": "test_disease_term_id"
                        }
                    ],
                    "donor_id": [
                        "test_donor_1"
                    ],
                    "id": mock.ANY,
                    "is_primary_data": "BOTH",
                    "is_valid": True,
                    "mean_genes_per_cell": 0.5,
                    "name": "test_dataset_name",
                    "organism": [
                        {
                            "label": "test_organism_label",
                            "ontology_term_id": "test_organism_term_id"
                        }
                    ],
                    "processing_status": {
                        "created_at": 0,
                        "cxg_status": "NA",
                        "dataset_id": mock.ANY,
                        "h5ad_status": "NA",
                        "id": "NA",
                        "processing_status": "PENDING",
                        "rds_status": "NA",
                        "updated_at": 0,
                        "upload_progress": 1,
                        "upload_status": "WAITING",
                        "validation_status": "NA"
                    },
                    "published": True,
                    "published_at": 1234,
                    "revision": 0, # NA
                    "schema_version": "3.0.0",
                    "self_reported_ethnicity": [
                        {
                            "label": "test_self_reported_ethnicity_label",
                            "ontology_term_id": "test_self_reported_ethnicity_term_id"
                        }
                    ],
                    "sex": [
                        {
                            "label": "test_sex_label",
                            "ontology_term_id": "test_sex_term_id"
                        }
                    ],
                    "suspension_type": [
                        "test_suspension_type"
                    ],
                    "tissue": [
                        {
                            "label": "test_tissue_label",
                            "ontology_term_id": "test_tissue_term_id"
                        }
                    ],
                    "tombstone": False,
                    "updated_at": 1234,
                    "x_approximate_distribution": "normal"
                },
                {
                    "assay": [
                        {
                            "label": "test_assay_label",
                            "ontology_term_id": "test_assay_term_id"
                        }
                    ],
                    "batch_condition": [
                        "test_batch_1",
                        "test_batch_2"
                    ],
                    "cell_count": 10,
                    "cell_type": [
                        {
                            "label": "test_cell_type_label",
                            "ontology_term_id": "test_cell_type_term_id"
                        }
                    ],
                    "collection_id": collection.collection_id.id,
                    "created_at": 1234,
                    "dataset_assets": [],
                    "dataset_deployments": [
                        {
                            "url": "TODO"
                        }
                    ],
                    "development_stage": [
                        {
                            "label": "test_development_stage_label",
                            "ontology_term_id": "test_development_stage_term_id"
                        }
                    ],
                    "disease": [
                        {
                            "label": "test_disease_label",
                            "ontology_term_id": "test_disease_term_id"
                        }
                    ],
                    "donor_id": [
                        "test_donor_1"
                    ],
                    "id": mock.ANY,
                    "is_primary_data": "BOTH",
                    "is_valid": True,
                    "mean_genes_per_cell": 0.5,
                    "name": "test_dataset_name",
                    "organism": [
                        {
                            "label": "test_organism_label",
                            "ontology_term_id": "test_organism_term_id"
                        }
                    ],
                    "processing_status": {
                        "created_at": 0,
                        "cxg_status": "NA",
                        "dataset_id": mock.ANY,
                        "h5ad_status": "NA",
                        "id": "NA",
                        "processing_status": "PENDING",
                        "rds_status": "NA",
                        "updated_at": 0,
                        "upload_progress": 1,
                        "upload_status": "WAITING",
                        "validation_status": "NA"
                    },
                    "published": True,
                    "published_at": 1234,
                    "revision": 0,
                    "schema_version": "3.0.0",
                    "self_reported_ethnicity": [
                        {
                            "label": "test_self_reported_ethnicity_label",
                            "ontology_term_id": "test_self_reported_ethnicity_term_id"
                        }
                    ],
                    "sex": [
                        {
                            "label": "test_sex_label",
                            "ontology_term_id": "test_sex_term_id"
                        }
                    ],
                    "suspension_type": [
                        "test_suspension_type"
                    ],
                    "tissue": [
                        {
                            "label": "test_tissue_label",
                            "ontology_term_id": "test_tissue_term_id"
                        }
                    ],
                    "tombstone": False,
                    "updated_at": 1234,
                    "x_approximate_distribution": "normal"
                }
            ],
            "description": "described",
            "id": mock.ANY,
            "links": [],
            "name": "test_collection",
            "published_at": 1234,
            "updated_at": 1234,
            "visibility": "PUBLIC"
        }

        with self.subTest("auth cookie"):
            expected_body["access_type"] = "WRITE"
            test_url = furl(path=f"/dp/v1/collections/{collection.version_id}", query_params=dict(visibility="PUBLIC"))
            cxguser_cookie = self.get_cxguser_token()
            response = self.app.get(test_url.url, headers=dict(host="localhost", Cookie=cxguser_cookie))
            self.assertEqual(200, response.status_code)
            actual_body = self.remove_timestamps(json.loads(response.data))
            self.assertDictEqual(actual_body, expected_body)

        with self.subTest("no auth cookie"):
            expected_body["access_type"] = "READ"
            test_url = furl(path=f"/dp/v1/collections/{collection.version_id}", query_params=dict(visibility="PUBLIC"))
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            self.assertEqual(200, response.status_code)
            actual_body = self.remove_timestamps(json.loads(response.data))
            self.assertDictEqual(actual_body, expected_body)

    # ✅
    def test__get_collection__ok(self):
        # Generate test cases
        authenticated = [True, False]
        owns = [True, False]
        visibility = ["public", "private"]
        test_cases = [params for params in itertools.product(authenticated, owns, visibility)]

        # Generate test collection
        # Note: for private collections, you want to use version_id
        test_collections = dict(
            public_not_owner=self.generate_published_collection(owner="someone else").collection_id.id,
            public_owned=self.generate_published_collection(owner="test_user_id").collection_id.id,
            private_not_owner=self.generate_unpublished_collection(owner="someone else").version_id.id,
            private_owned=self.generate_unpublished_collection(owner="test_user_id").version_id.id,
        )

        # run
        for auth, owns, visi in test_cases:
            expected_response_code = 200

            test_collection_id = test_collections["_".join([visi, "owned" if owns else "not_owner"])]
            expected_access_type = "WRITE" if owns and auth else "READ"

            with self.subTest(f"auth:{auth}, owns:{owns}, visi:{visi}, acc:{expected_access_type}"):
                test_url = furl(path=f"/dp/v1/collections/{test_collection_id}")

                headers = dict(host="localhost")
                if auth:
                    headers["Cookie"] = self.get_cxguser_token()
                print("headers - ", headers)
                response = self.app.get(test_url.url, headers=headers)

                self.assertEqual(expected_response_code, response.status_code)
                if expected_response_code == 200:
                    actual_body = json.loads(response.data)
                    self.assertEqual(expected_access_type, actual_body["access_type"])

    # ✅
    def test__get_collection_id__403_not_found(self):
        """Verify the test collection exists and the expected fields exist."""
        test_url = furl(path="/dp/v1/collections/AAAA-BBBB-CCCC-DDDD", query_params=dict(visibility="PUBLIC"))
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, response.status_code)

    # ✅, but this test should do more assertions
    def test__post_collection_returns_id_on_success(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)

        # Add curator_name
        data["curator_name"] = "john smith"
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)

    # ✅
    def test__post_collection_normalizes_doi(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [
                {"link_name": "DOI Link", "link_url": "10.1016/foo", "link_type": "DOI"},
            ],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        collection = self.business_logic.get_collection_version(CollectionVersionId(collection_id))
        print(collection)
        doi = next(link.uri for link in collection.metadata.links if link.type == "DOI") # TODO: careful
        self.assertEquals(doi, "https://doi.org/10.1016/foo")

    # ✅
    def test__post_collection_rejects_two_dois(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
            ],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(400, response.status_code)

    # ✅
    def test__post_collection_rejects_doi_not_in_crossref(self):
        self.crossref_provider.fetch_metadata = Mock(
            side_effect=CrossrefDOINotFoundException("Mocked CrossrefDOINotFoundException")
        )
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(400, response.status_code)
        error_payload = json.loads(response.data)
        self.assertEqual(error_payload["detail"][0], {"link_type": "DOI", "reason": "DOI cannot be found on Crossref"})


    # ✅
    def test__post_collection_rejects_invalid_doi(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [{"link_name": "DOI Link", "link_url": "invalid/doi", "link_type": "DOI"}],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(400, response.status_code)
        error_payload = json.loads(response.data)
        self.assertEqual([{"link_type": "DOI", "reason": "Invalid DOI"}], error_payload["detail"])

    # ✅
    def test__post_collection_ignores_metadata_if_crossref_exception(self):
        self.crossref_provider.fetch_metadata = Mock(
            side_effect=CrossrefFetchException("Mocked CrossrefFetchException")
        )
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        collection = self.business_logic.get_collection_version(CollectionVersionId(collection_id))
        self.assertIsNone(collection.publisher_metadata)

    # ✅
    def test__post_collection_ignores_metadata_if_no_doi(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        collection = self.business_logic.get_collection_version(CollectionVersionId(collection_id))
        self.assertIsNone(collection.publisher_metadata)

    # ✅
    def test__post_collection_adds_publisher_metadata(self):

        self.crossref_provider.fetch_metadata = Mock(
            return_value=generate_mock_publisher_metadata()
        )

        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        collection = self.business_logic.get_collection_version(CollectionVersionId(collection_id))
        self.assertIsNotNone(collection.publisher_metadata)
        self.assertEqual(collection.publisher_metadata, generate_mock_publisher_metadata())

    # ✅
    def test__post_collection_fails_with_extra_fields(self):
        test_url = furl(path="/dp/v1/collections")
        test_data = [
            {
                "name": "extra field in collection",
                "owner": "someone else",
                "description": "This is a test collection",
                "contact_name": "person human",
                "curator_name": "",
                "contact_email": "person@human.com",
                "links": [
                    {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                    {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
                ],
            },
            {
                "name": "extra field in link",
                "description": "This is a test collection",
                "contact_name": "person human",
                "curator_name": "",
                "contact_email": "person@human.com",
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
                        "Cookie": self.get_cxguser_token(),
                    },
                    data=json.dumps(data),
                )
                self.assertEqual(400, response.status_code)

    # ✅
    def test__post_collection_fails_if_data_missing(self):
        test_url = furl(path="/dp/v1/collections")
        data = json.dumps({"name": "bkjbjbjmbjm"})
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=data,
        )
        self.assertEqual(400, response.status_code)

    # ✅
    def test__can_retrieve_created_collection(self):
        test_url = furl(path="/dp/v1/collections")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        data = {
            "name": "another collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "curator_name": "",
            "contact_email": "person@human.com",
            "links": [
                {"link_url": "http://doi.org/10.1016", "link_type": "OTHER"},
                {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
            ],
        }
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(data))
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]

        test_url = furl(path=f"/dp/v1/collections/{collection_id}")
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
        test_url = furl(path=f"/dp/v1/collections/{collection_id}")
        response = self.app.get(test_url.url, headers=no_cookie_headers)
        self.assertEqual("READ", json.loads(response.data)["access_type"])

    # 🔴 used to be done by the database - we need to decide where this goes now
    def test__create_collection_strip_string_fields(self):
        test_url = furl(path="/dp/v1/collections")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        data = {
            "name": "   another collection name   ",
            "description": "    This is a test collection  ",
            "contact_name": " person human   ",
            "contact_email": "person@human.com  ",
            "curator_name": "",
            "links": [
                {"link_url": "     http://doi.org/10.1016  ", "link_type": "OTHER"},
                {"link_name": "  DOI Link 2", "link_url": "http://doi.org/10.1017   ", "link_type": "DOI"},
            ],
        }
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(data))
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]

        test_url = furl(path=f"/dp/v1/collections/{collection_id}")
        test_url.add(query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        self.assertEqual(body["description"], data["description"].strip())
        self.assertEqual(body["name"], data["name"].strip())
        self.assertEqual(body["contact_name"], body["contact_name"].strip())
        self.assertEqual(body["contact_email"], body["contact_email"].strip())
        self.assertEqual(body["data_submission_policy_version"], body["data_submission_policy_version"].strip())

        for link in body["links"]:
            self.assertEqual(link["link_url"], link["link_url"].strip())

    # 🔴 TODO: this needs attention
    def test__list_collection__check_owner(self):

        # Generate test collection
        public_owned = self.generate_published_collection(owner="test_user_id").collection_id
        private_owned = self.generate_unpublished_collection(owner="test_user_id").version_id
        public_not_owned = self.generate_published_collection(owner="someone else").collection_id
        private_not_owned = self.generate_unpublished_collection(owner="someone else").version_id

        revision_not_owned = self.business_logic.create_collection_version(public_not_owned).version_id
        revision_owned = self.business_logic.create_collection_version(public_owned).version_id

        path = "/dp/v1/collections"
        with self.subTest("no auth"):
            headers = {"host": "localhost", "Content-Type": "application/json"}
            response = self.app.get(path, headers=headers)
            self.assertEqual(200, response.status_code)
            result = json.loads(response.data)
            collections = result.get("collections")
            self.assertIsNotNone(collections)
            ids = [collection.get("id") for collection in collections]
            self.assertIn(public_owned.id, ids)
            self.assertIn(public_not_owned.id, ids)
            self.assertNotIn(private_owned.id, ids)
            self.assertNotIn(private_not_owned.id, ids)
            self.assertNotIn(revision_owned.id, ids)

        with self.subTest("auth"):
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
            response = self.app.get(path, headers=headers)
            self.assertEqual(200, response.status_code)
            result = json.loads(response.data)
            collections = result.get("collections")
            self.assertIsNotNone(collections)
            ids = [collection.get("id") for collection in collections]
            self.assertIn(public_owned.id, ids)
            self.assertIn(public_not_owned.id, ids)
            self.assertIn(private_owned.id, ids)
            self.assertNotIn(private_not_owned.id, ids)
            self.assertNotIn(revision_not_owned.id, ids)
            self.assertIn(revision_owned.id, ids)
            self.assertTrue(
                [collection for collection in collections if collection.get("revision_of") == public_owned][0]
            )

    # ✅
    def test__get_all_collections_for_index(self):
        """
        The `collections/index` endpoint should only return public collections
        """

        collection = self.generate_published_collection()
        private_collection = self.generate_unpublished_collection()

        # TODO: a tombstoned collection should not be returned as well

        test_url = furl(path="/dp/v1/collections/index")
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [collection["id"] for collection in body]
        self.assertIn(collection.collection_id.id, ids)
        self.assertNotIn(private_collection.collection_id.id, ids)
        self.assertNotIn(private_collection.version_id.id, ids)

        actual_collection = body[-1]  # last added collection
        self.assertEqual(actual_collection["id"], collection.collection_id.id)
        self.assertEqual(actual_collection["name"], collection.metadata.name)
        self.assertNotIn("description", actual_collection)
        # Both `published_at` and `revised_at` should point to the same timestamp
        self.assertEqual(actual_collection["published_at"], collection.published_at.timestamp())
        self.assertEqual(actual_collection["revised_at"], collection.published_at.timestamp())

    # ✅
    def test__create_collection__InvalidParameters_DOI(self):
        tests = [
            (
                dict(
                    name="not blank",
                    description="description",
                    contact_name="some name",
                    contact_email="robot@email.com",
                    links=[{"link_type": "DOI", "link_url": "bad_doi"}],
                ),
                [
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
                headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
                response = self.app.post("/dp/v1/collections", headers=headers, data=json.dumps(body))
                self.assertEqual(400, response.status_code)
                for error in expected_errors:
                    self.assertIn(error, response.json["detail"])


# 🔴 TODO: This should be reviewed. Collection deletion is a weird beast
class TestCollectionDeletion(NewBaseTest):
    def test_delete_private_collection__ok(self):
        # Generate test collection
        collection = self.generate_unpublished_collection()
        processing_status_1 = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        processing_status_2 = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 100.0}

        dataset_1 = self.generate_dataset_with_s3_resources(
            self.session, collection=collection, processing_status=processing_status_1
        )
        dataset_2 = self.generate_dataset_with_s3_resources(
            self.session, collection=collection, processing_status=processing_status_2
        )

        s3_objects = self.get_s3_object_paths_from_dataset(dataset_1) + self.get_s3_object_paths_from_dataset(dataset_2)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
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

    def test_delete_collection_revision__ok(self):
        # Generate test collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )
        # Generate the public collection with the same id as the private so a tombstone is created
        revision = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id", revision_of=collection.id
        )

        processing_status_1 = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        processing_status_2 = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 100.0}

        dataset_1 = self.generate_dataset(self.session, collection=revision, processing_status=processing_status_1)
        dataset_2 = self.generate_dataset(self.session, collection=revision, processing_status=processing_status_2)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        test_private_url = furl(path=f"/dp/v1/collections/{revision.id}")
        test_public_url = furl(path=f"/dp/v1/collections/{collection.id}")
        response = self.app.get(test_private_url.url, headers=headers)
        self.assertEqual(200, response.status_code)

        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset_1.id, dataset_ids)
        self.assertIn(dataset_2.id, dataset_ids)

        # delete collection
        response = self.app.delete(test_private_url.url, headers=headers)

        self.assertEqual(response.status_code, 204)

        # check collection and datasets delete
        response = self.app.get(test_private_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

        # Public collection still exists
        response = self.app.get(test_public_url.url, headers=headers)
        self.assertEqual(response.status_code, 200)

    def test_tombstone_published_collection_with_revision__ok(self):
        """Both the published and revised collections are tombstoned."""
        # Generate the public collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )
        # Generate test collection
        revision = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id", revision_of=collection.id
        )
        revision_id = revision.id

        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 100.0}

        dataset_rev = self.generate_dataset(self.session, collection=revision, processing_status=processing_status)
        dataset_pub = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}

        # Verify private collections exist
        test_private_url = furl(path=f"/dp/v1/collections/{revision.id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_private_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset_rev.id, dataset_ids)

        # Verify public collections exist
        test_public_url = furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PUBLIC"))
        response = self.app.get(test_public_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset_pub.id, dataset_ids)

        # delete public collection
        response = self.app.delete(test_public_url.url, headers=headers)
        self.assertEqual(response.status_code, 204)

        # check collection revision and datasets are gone
        response = self.app.get(test_private_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

        # check public collection is tombstoned and datasets deleted.
        response = self.app.get(test_public_url.url, headers=headers)
        self.assertEqual(response.status_code, 410)

        self.session.expire_all()
        collection = Collection.get_collection(
            self.session, collection.id, CollectionVisibility.PUBLIC.name, include_tombstones=True
        )
        self.assertTrue(collection.tombstone)
        self.assertTrue(dataset_pub.tombstone)
        rev_collection = Collection.get_collection(self.session, revision_id, include_tombstones=True)
        self.assertIsNone(rev_collection)  # Revision should be deleted, not tombstoned

    def test_delete_collection__dataset_not_available(self):
        # Generate the public collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id", tombstone=True
        )
        # Generate test collection
        revision = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id", revision_of=collection.id
        )
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 100.0}

        dataset = self.generate_dataset(self.session, collection=revision, processing_status=processing_status)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        dataset_url = furl(path=f"/dp/v1/datasets/{dataset.id}/status")
        response = self.app.get(dataset_url.url, headers=headers)
        self.assertEqual(response.status_code, 200)

        test_url = furl(path=f"/dp/v1/collections/{revision.id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.delete(test_url.url, headers=headers)

        self.assertEqual(response.status_code, 204)

        response = self.app.get(dataset_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_delete_collection__already_tombstoned__ok(self):
        # Generate the public collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id", tombstone=True
        )
        # Generate test collection
        self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id", revision_of=collection.id
        )

        test_url = furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PRIVATE"))
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_delete_collection__public__ok(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )

        test_urls = [
            furl(path=f"/dp/v1/collections/{collection.id}"),
        ]
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        for test_url in test_urls:
            with self.subTest(test_url.url):
                response = self.app.delete(test_url.url, headers=headers)
                self.assertEqual(response.status_code, 204)

                response = self.app.get(test_url.url, headers=headers)
                body = json.loads(response.data)
                self.assertEqual(response.status_code, 410)
                self.assertEqual(body, "")

    def test_delete_collection__not_owner(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="someone_else"
        )
        test_url = furl(path=f"/dp/v1/collections/{collection.id}")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_delete_collection__does_not_exist(self):
        fake_id = generate_id()
        test_url = furl(path=f"/dp/v1/collections/{fake_id}", query_params=dict(visibility="PRIVATE"))
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
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
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.get("/dp/v1/collections", headers=headers)

        print(f"%%%%%%%%%%%% {response.data}")

        collection_ids = [collection["id"] for collection in json.loads(response.data)["collections"]]
        self.assertIn(private_collection.id, collection_ids)
        self.assertIn(public_collection.id, collection_ids)
        self.assertIn(collection_to_delete.id, collection_ids)

        test_url = furl(path=f"/dp/v1/collections/{collection_to_delete.id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 204)

        # check not returned privately
        response = self.app.get("/dp/v1/collections", headers=headers)
        collection_ids = [collection["id"] for collection in json.loads(response.data)["collections"]]
        self.assertIn(private_collection.id, collection_ids)
        self.assertIn(public_collection.id, collection_ids)

        self.assertNotIn(collection_to_delete.id, collection_ids)

        # check not returned publicly
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get("/dp/v1/collections", headers=headers)
        collection_ids = [collection["id"] for collection in json.loads(response.data)["collections"]]
        self.assertIn(public_collection.id, collection_ids)
        self.assertNotIn(private_collection.id, collection_ids)
        self.assertNotIn(collection_to_delete.id, collection_ids)


class TestUpdateCollection(NewBaseTest):

    # ✅
    def test__update_collection__OK(self):
        collection = self.generate_unpublished_collection()
        test_fields = [
            "name",
            "description",
            "contact_name",
            "contact_email",
            "links",
        ]
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # Update the collection
        expected_body = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
        }
        data = json.dumps(expected_body)
        response = self.app.put(f"/dp/v1/collections/{collection.version_id.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        for field in test_fields:
            self.assertEqual(expected_body[field], actual_body[field])

    # ✅
    def test__update_collection__403(self):
        collection = self.generate_unpublished_collection(owner="someone else")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        data = json.dumps({"name": "new name"})
        response = self.app.put(f"/dp/v1/collections/{collection.version_id.id}", data=data, headers=headers)
        self.assertEqual(403, response.status_code)

    # ✅
    def test__update_collection_links__OK(self):
        collection = self.generate_unpublished_collection()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # add links
        links = [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}]
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.version_id.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])

        # remove links
        links.pop()
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.version_id.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])

        # update links
        links = [{"link_name": "New name", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}]
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.version_id.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])

        # all together
        links = [
            {"link_name": "Link 1", "link_url": "http://link.com", "link_type": "OTHER"},
            {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
        ]
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.version_id.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])

        # Clear All links
        links = []
        data = json.dumps({"links": links})
        response = self.app.put(f"/dp/v1/collections/{collection.version_id}", data=data, headers=headers)

        self.assertEqual(200, response.status_code)
        self.assertEqual(links, json.loads(response.data)["links"])

    # ✅
    def test__update_collection_new_doi_updates_metadata(self):

        # Generate a collection with "Old Journal" as publisher metadata
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata("Old Journal"))
        collection = self.generate_unpublished_collection(
            links=[Link("Link 1", "DOI", "http://doi.org/123")]
        )

        self.assertIsNotNone(collection.publisher_metadata)
        if collection.publisher_metadata: # pylance
            self.assertEqual("Old Journal", collection.publisher_metadata["journal"]) 
        
        # From now on, Crossref will return `New Journal`
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata("New Journal"))

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id}",
            data=json.dumps({"links": [{"link_name": "Link 1", "link_url": "http://doi.org/456", "link_type": "DOI"}]}),
            headers=headers,
        )
        self.assertEqual(200, response.status_code)

        self.crossref_provider.fetch_metadata.assert_called_once()

        actual_body = json.loads(response.data)
        self.assertIsNotNone(actual_body["publisher_metadata"])
        self.assertIsNotNone(actual_body["publisher_metadata"]["journal"])
        self.assertEqual("New Journal", actual_body["publisher_metadata"]["journal"])

    # ✅
    def test__update_collection_remove_doi_deletes_metadata(self):

        # Generate a collection with "Old Journal" as publisher metadata
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata("Old Journal"))
        collection = self.generate_unpublished_collection(
            links=[Link("Link 1", "DOI", "http://doi.org/123")]
        )

        self.assertIsNotNone(collection.publisher_metadata)
        if collection.publisher_metadata: # pylance
            self.assertEqual("Old Journal", collection.publisher_metadata["journal"]) 

        # From now on, Crossref will return `New Journal`
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata("New Journal"))

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        # We're passing an empty links object, therefore the DOI is deleted
        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id}",
            data=json.dumps({"links": []}),
            headers=headers,
        )

        self.assertEqual(200, response.status_code)
        # No Crossref calls should happen
        self.crossref_provider.fetch_metadata.assert_not_called()

        # The `publisher_metadata` node should not longer be in the collection
        actual_body = json.loads(response.data)
        self.assertNotIn("publisher_metadata", actual_body)

    # ✅
    def test__update_collection_same_doi_does_not_update_metadata(self):
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata("Old Journal"))

        collection = self.generate_unpublished_collection(
            links=[Link("Link 1", "DOI", "http://doi.org/123")]
        )

        self.assertIsNotNone(collection.publisher_metadata)
        if collection.publisher_metadata: # pylance
            self.assertEqual("Old Journal", collection.publisher_metadata["journal"]) 

        # From now on, Crossref will return `New Journal`
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata("New Journal"))

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        # Note that the DOI is the same as the original
        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id}",
            data=json.dumps({"links": [{"link_name": "Link 1", "link_url": "http://doi.org/123", "link_type": "DOI"}]}),
            headers=headers,
        )

        self.assertEqual(200, response.status_code)
        # No Crossref calls should happen, since the DOI is unmodified
        self.crossref_provider.fetch_metadata.assert_not_called()

        # The `publisher_metadata` node should exist and be the same
        actual_body = json.loads(response.data)
        self.assertIsNotNone(actual_body["publisher_metadata"])
        self.assertIsNotNone(actual_body["publisher_metadata"]["journal"])
        self.assertEqual("Old Journal", actual_body["publisher_metadata"]["journal"])


class TestCollectionsCurators(NewBaseTest):


    # 🔴 not quite sure why this fails - a non curator should not have access_type of WRITE
    def test_view_non_owned_private_collection__ok(self):
        # Generate test collection
        collection = self.generate_unpublished_collection(owner="another_test_user_id")
        

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        test_url = furl(path=f"/dp/v1/collections/{collection.version_id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers=headers)

        # This will pass even for non curators.
        # Why are users allowed to view private collections that they don't own?
        self.assertEqual(response.status_code, 200)

        body = json.loads(response.data)
        self.assertEqual(body["access_type"], "WRITE")

    # ✅
    def test_update_non_owned_private_collection_as_super_curator__ok(self):
        # Generate test collection
        collection = self.generate_unpublished_collection(owner="another_test_user_id")

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("super")}

        modified_collection = {
            "name": "new name",
            "description": collection.metadata.description,
            "contact_name": collection.metadata.contact_name,
            "contact_email": collection.metadata.contact_email,
            "links": collection.metadata.links,
        }
        data = json.dumps(modified_collection)
        response = self.app.put(f"/dp/v1/collections/{collection.version_id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)

    # 🔴 - delete collection will be implemented later
    def test_delete_non_owned_private_collection_as_super_curator__ok(self):
        # Generate test collection
        collection = self.generate_unpublished_collection(owner="another_test_user_id")

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("super")}

        test_url = furl(path=f"/dp/v1/collections/{collection.version_id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(204, response.status_code)


# TODO: 💛 Not an API test, but still valuable. Figure out where to put it. Maybe create a test_validation.py?
class TestVerifyCollection(unittest.TestCase):
    def test_empty_body(self):
        body = dict()
        errors = []
        verify_collection_body(body, errors)
        self.assertFalse(errors)

    def test_blank_fields(self):
        errors = []
        body = dict(name="", contact_name="", description="", contact_email="")
        verify_collection_body(body, errors)
        error_message = "Cannot be blank."
        self.assertIn({"name": "description", "reason": error_message}, errors)
        self.assertIn({"name": "name", "reason": error_message}, errors)
        self.assertIn({"name": "contact_name", "reason": error_message}, errors)
        self.assertIn({"name": "contact_email", "reason": error_message}, errors)

    def test_invalid_characters_in_field(self):
        invalid_strings = [b"\x00some data", b"text\x1f", b"text\x01", b"\x7ftext"]
        for test_string in invalid_strings:
            with self.subTest(test_string):
                errors = []
                string = test_string.decode(encoding="utf-8")
                body = dict(name=string, contact_name=string, description=string, contact_email="email@email.com")
                verify_collection_body(body, errors)
                error_message = "Invalid characters detected."
                self.assertEqual(1, len(errors))
                self.assertIn({"name": "name", "reason": error_message}, errors)

    def test_invalid_email(self):
        bad_emails = ["@.", "email@.", "@place.com", "email@.com", "email@place."]
        body = dict()

        for email in bad_emails:
            with self.subTest(email):
                body["contact_email"] = email
                errors = []
                verify_collection_body(body, errors)
                self.assertEqual([{"name": "contact_email", "reason": "Invalid format."}], errors)

    def test_OK(self):
        body = dict(name="something", contact_name="a name", description="description", contact_email="email@place.com")
        errors = []
        verify_collection_body(body, errors)
        self.assertFalse(errors)

    def test__link__INVALID(self):
        test_urls = ["://", "google", ".com", "google.com", "https://"]
        for link_type in ProjectLinkType:
            if link_type.name == "DOI":
                continue
            for test_url in test_urls:
                link_body = [{"link_type": link_type.name, "link_url": test_url}]
                with self.subTest(link_body):
                    errors = []
                    body = dict(links=link_body)
                    verify_collection_body(body, errors)
                    expected_error = [dict(reason="Invalid URL.", name="links[0]", value=link_body[0]["link_url"])]
                    self.assertEqual(expected_error, errors)

    def test__link__OK(self):
        test_urls = ["https://www.google.com", "http://somewhere.org/path/?abcd=123"]
        for link_type in ProjectLinkType:
            if link_type.name == "DOI":
                continue
            for test_url in test_urls:
                link_body = [{"link_type": link_type.name, "link_url": test_url}]
                with self.subTest(link_body):
                    errors = []
                    body = dict(links=link_body)
                    verify_collection_body(body, errors)
                    self.assertFalse(errors)


# TODO: these tests all require the generation of a dataset
class TestDataset(NewBaseTest):

    # ✅
    def test__post_dataset_asset__OK(self):
        self.business_logic.get_dataset_artifact_download_data = Mock(
            return_value=DatasetArtifactDownloadData("asset.h5ad", "H5AD", 1000, "http://presigned.url")
        )
        version = self.generate_dataset(artifacts=[DatasetArtifactUpdate("H5AD", "http://mock.uri/asset.h5ad")])
        dataset_version_id = version.dataset_version_id
        artifact_id = version.artifact_ids[0]

        test_url = furl(path=f"/dp/v1/datasets/{dataset_version_id}/asset/{artifact_id}")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(200, response.status_code)

        actual_body = json.loads(response.data)
        self.assertEqual(actual_body["presigned_url"], "http://presigned.url")
        self.assertEqual(actual_body["dataset_id"], version.dataset_version_id)
        self.assertEqual(actual_body["file_size"], 1000)
        self.assertEqual(actual_body["file_name"], "asset.h5ad")


    # ✅, but I think the behavior could be improved
    def test__post_dataset_asset__file_SERVER_ERROR(self):
        """
        `post_dataset_asset` should throw 500 if presigned_url or file_size aren't returned from the server
        """
        version = self.generate_dataset(artifacts=[DatasetArtifactUpdate("H5AD", "http://mock.uri/asset.h5ad")])
        dataset_version_id = version.dataset_version_id
        artifact_id = version.artifact_ids[0]

        test_url = furl(path=f"/dp/v1/datasets/{dataset_version_id}/asset/{artifact_id}")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(500, response.status_code)
        body = json.loads(response.data)
        self.assertEqual("An internal server error has occurred. Please try again later.", body["detail"])

    # ✅
    def test__post_dataset_asset__dataset_NOT_FOUND(self):
        test_url = furl(path="/dp/v1/datasets/test_user_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.data)
        self.assertEqual("'dataset/test_user_id' not found.", body["detail"])

    # ✅
    def test__post_dataset_asset__asset_NOT_FOUND(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/fake_asset")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.data)
        self.assertEqual("'dataset/test_dataset_id/asset/fake_asset' not found.", body["detail"])

    # ✅
    def test__get_status__ok(self):
        dataset = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate("processing_status", DatasetProcessingStatus.PENDING),
                DatasetStatusUpdate("upload_status", DatasetUploadStatus.UPLOADING),
            ]
        )
        # TODO: why do we need processing_status_id? we can probably remove
        test_url = furl(path=f"/dp/v1/datasets/{dataset.dataset_version_id}/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        expected_body = {
            "cxg_status": "NA",
            "rds_status": "NA",
            "h5ad_status": "NA",
            "processing_status": "PENDING",
            "dataset_id": dataset.dataset_version_id,
            "id": "NA", # TODO: I am deprecating this, I don't think it has any use.
            "upload_progress": 1,
            "upload_status": "UPLOADING",
            "validation_status": "NA",
        }
        self.assertEqual(expected_body, actual_body)

    # ✅
    def test__get_status__403(self):
        dataset = self.generate_dataset(
            owner="someone_else",
        )
        test_url = furl(path=f"/dp/v1/datasets/{dataset.dataset_version_id}/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(403, response.status_code)

    # 💛, passes, but review the assertions
    def test__get_all_datasets_for_index(self):
        # TODO: here we probably don't need to test the logic - needs to simplify this

        private_dataset = self.generate_dataset()
        public_dataset = self.generate_dataset(publish=True)

        test_url = furl(path="/dp/v1/datasets/index")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [d["id"] for d in body]
        self.assertIn(public_dataset.dataset_id, ids)
        self.assertNotIn(private_dataset.dataset_id, ids)

        actual_dataset = None
        for d in body:
            if d["id"] == public_dataset.dataset_id:
                actual_dataset = d

        persisted_dataset = self.business_logic.get_dataset_version(DatasetVersionId(public_dataset.dataset_version_id))
        self.assertIsNotNone(actual_dataset)
        self.assertIsNotNone(persisted_dataset)

        def convert_ontology(ontologies):
            return [dataclasses.asdict(o) for o in ontologies]

        if actual_dataset is not None and persisted_dataset is not None: #pylance

            self.assertNotIn("description", actual_dataset)
            self.assertEqual(actual_dataset["id"], persisted_dataset.dataset_id.id)
            self.assertEqual(actual_dataset["name"], persisted_dataset.metadata.name)
            # self.assertEqual(actual_dataset["collection_id"], persisted_dataset.collection_id)
            self.assertEqual(actual_dataset["assay"], convert_ontology(persisted_dataset.metadata.assay))
            self.assertEqual(actual_dataset["tissue"], convert_ontology(persisted_dataset.metadata.tissue))
            self.assertEqual(actual_dataset["disease"], convert_ontology(persisted_dataset.metadata.disease))
            self.assertEqual(actual_dataset["sex"], convert_ontology(persisted_dataset.metadata.sex))
            self.assertEqual(actual_dataset["self_reported_ethnicity"], convert_ontology(persisted_dataset.metadata.self_reported_ethnicity))
            self.assertEqual(actual_dataset["organism"], convert_ontology(persisted_dataset.metadata.organism))
            self.assertEqual(actual_dataset["development_stage"], convert_ontology(persisted_dataset.metadata.development_stage))
            self.assertEqual(actual_dataset["cell_count"], persisted_dataset.metadata.cell_count)
            self.assertEqual(actual_dataset["cell_type"], convert_ontology(persisted_dataset.metadata.cell_type))
            self.assertEqual(actual_dataset["is_primary_data"], persisted_dataset.metadata.is_primary_data)
            self.assertEqual(actual_dataset["mean_genes_per_cell"], persisted_dataset.metadata.mean_genes_per_cell)
            # self.assertEqual(actual_dataset["explorer_url"], persisted_dataset.explorer_url)
            # self.assertEqual(actual_dataset["published_at"], persisted_dataset.published_at.timestamp())
            # self.assertEqual(actual_dataset["revised_at"], persisted_dataset.revised_at.timestamp())

    # ✅
    def test__get_all_datasets_for_index_with_ontology_expansion(self):

        import copy
        modified_metadata = copy.deepcopy(self.sample_dataset_metadata)
        modified_metadata.development_stage = [OntologyTermId("Test", "HsapDv:0000008")]
        modified_metadata.tissue = [OntologyTermId("Test", "UBERON:0002048")]
        modified_metadata.cell_type = [OntologyTermId("Test", "CL:0000738")]

        dataset = self.generate_dataset(metadata=modified_metadata, publish=True)

        test_url = furl(path="/dp/v1/datasets/index")

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        actual_dataset = None
        for d in body:
            if d["id"] == dataset.dataset_id:
                actual_dataset = d
        self.assertIsNotNone(actual_dataset)

        def convert_ontology(ontologies):
            return [dataclasses.asdict(o) for o in ontologies]

        if actual_dataset is not None: # pylance

            self.assertEqual(actual_dataset["development_stage"], convert_ontology(modified_metadata.development_stage))
            self.assertEqual(
                actual_dataset["development_stage_ancestors"],
                ["HsapDv:0000008", "HsapDv:0000006", "HsapDv:0000002", "HsapDv:0000045", "HsapDv:0000001"],
            )

            self.assertEqual(actual_dataset["tissue"], convert_ontology(modified_metadata.tissue))
            self.assertCountEqual(
                actual_dataset["tissue_ancestors"],
                [
                    "UBERON:0001004",
                    "UBERON:0001005",
                    "UBERON:0000065",
                    "UBERON:0000170",
                    "UBERON:0002048",
                    "UBERON:0001558",
                    "UBERON:0000072",
                    "UBERON:0000171",
                ],
            )

            self.assertEqual(actual_dataset["cell_type"], convert_ontology(modified_metadata.cell_type))
            self.assertCountEqual(
                actual_dataset["cell_type_ancestors"],
                [
                    "CL:0000255",
                    "CL:0002371",
                    "CL:0000988",
                    "CL:0000738",
                    "CL:0000548",
                    "CL:0000219",
                    "CL:0000003",
                    "CL:0002242",
                ],
            )

    # ✅
    def test__get_dataset_assets(self):
        # TODO: I don't think `filename` is relevant - review
        dataset = self.generate_dataset(artifacts=[
            DatasetArtifactUpdate("CXG", "s3://mock-bucket/mock-key.cxg"),
            DatasetArtifactUpdate("H5AD", "s3://mock-bucket/mock-key.h5ad"),
        ])

        test_url = furl(path=f"/dp/v1/datasets/{dataset.dataset_version_id}/assets")
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        self.assertIn("assets", body)
        assets = body["assets"]
        self.assertEqual(len(assets), 2)
        self.assertEqual(assets[0]["s3_uri"], "s3://mock-bucket/mock-key.cxg")
        self.assertEqual(assets[1]["s3_uri"], "s3://mock-bucket/mock-key.h5ad")

    # ✅
    def test__cancel_dataset_download__ok(self):
        # Test pre upload
        # TODO: this might need additional business logic
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("upload_status", DatasetUploadStatus.WAITING)]
        )
        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # Test while uploading
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("upload_status", DatasetUploadStatus.UPLOADING)]
        )
        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

    # ✅
    def test__cancel_dataset_download__dataset_does_not_exist(self):
        test_url = "/dp/v1/datasets/missing_dataset_id"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    # ✅
    def test__delete_uploaded_dataset__ok(self):
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("upload_status", DatasetUploadStatus.UPLOADING)]
        )
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # check dataset in collection
        collection_url = furl(path=f"/dp/v1/collections/{dataset.collection_version_id}")
        # TODO: do we still have a query for private/public?
        response = self.app.get(collection_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset.dataset_id, dataset_ids)

        # delete dataset
        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

        # check dataset no longer returned in collection
        response = self.app.get(collection_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertNotIn(dataset.dataset_version_id, dataset_ids)

    # ✅
    def test__call_delete_dataset__twice(self):
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("upload_status", DatasetUploadStatus.UPLOADING)]
        )
        # TODO: set upload progress to 10.0

        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # delete again
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    # ✅
    def test__delete_public_dataset_returns__405(self):

        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("upload_status", DatasetUploadStatus.UPLOADED)],
            publish=True,
        )

        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(405, response.status_code)
        self.assertEqual("Cannot delete a public Dataset", json.loads(response.data)["detail"])

    # ✅
    def test__cancel_dataset_download__user_not_collection_owner(self):

        dataset = self.generate_dataset(
            owner="someone_else",
            statuses=[DatasetStatusUpdate("upload_status", DatasetUploadStatus.WAITING)],
        )

        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    # ✅
    def test__cancel_dataset_download__user_not_logged_in(self):

        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("upload_status", DatasetUploadStatus.WAITING)],
        )

        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 401)

    # ✅
    def test__dataset_meta__ok(self):

        headers = {"host": "localhost", "Content-Type": "application/json"}

        with self.subTest("dataset is public"):

            test_uri_0 = "some_uri_0"
            
            public_dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate("CXG", test_uri_0)],
                publish=True,
            )

            test_url_public = f"/dp/v1/datasets/meta?url={public_dataset.explorer_url}"

            response = self.app.get(test_url_public, headers)
            self.assertEqual(response.status_code, 200)

            expected_identifiers = {
                "s3_uri": test_uri_0,
                "dataset_id": public_dataset.dataset_id,
                "collection_id": public_dataset.collection_id,
                "collection_visibility": "PUBLIC", # this is a published collection
                "tombstoned": False,
            }

            self.assertEqual(json.loads(response.data), expected_identifiers)

        with self.subTest("dataset is private"):
            test_uri_1 = "some_uri_1"
   
            private_dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate("CXG", test_uri_1)],
                publish=False,
            )

            test_url_private = f"/dp/v1/datasets/meta?url={private_dataset.explorer_url}"
            expected_identifiers = {
                "s3_uri": test_uri_1,
                "dataset_id": private_dataset.dataset_id,
                "collection_id": private_dataset.collection_id,
                "collection_visibility": "PRIVATE",
                "tombstoned": False,
            }

            response = self.app.get(test_url_private, headers)
            self.assertEqual(response.status_code, 200)

            self.assertEqual(json.loads(response.data), expected_identifiers)

    # ✅
    def test__dataset_meta__404(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}
        test_url_404 = "/dp/v1/datasets/meta?url=not_real"

        response = self.app.get(test_url_404, headers)
        self.assertEqual(response.status_code, 404)


class TestDatasetCurators(NewBaseTest):
    def setUp(self):
        super().setUp()

    def tearDown(self):
        super().tearDown()

    # 💛 review later
    def test__get_status__200_for_non_owned_dataset_as_super_curator(self):
        processing_status = DatasetStatusUpdate("upload_status", DatasetUploadStatus.WAITING)
        dataset = self.generate_dataset(owner="someone_else", statuses=[processing_status])
        test_url = f"/dp/v1/datasets/{dataset.dataset_id}/status"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("super")}
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(200, response.status_code)

    # 💛 review later
    def test__cancel_dataset_download__202_user_not_collection_owner_as_super_curator(self):
        processing_status = DatasetStatusUpdate("upload_status", DatasetUploadStatus.WAITING)
        dataset = self.generate_dataset(owner="someone_else", statuses=[processing_status])
        test_url = f"/dp/v1/datasets/{dataset.dataset_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("super")}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)



##### REVISIONS START HERE #####

class TestRevision(NewBaseTest):
    """Test case for starting a collection's revision."""

    # ✅
    def test__start_revision_of_a_collection__201(self):
        """
        Start a revision of a collection.
        """

        # Create a published collection with 2 datasets
        published_collection = self.generate_published_collection(add_datasets=2)

        # Starts a revision
        path = f"/dp/v1/collections/{published_collection.collection_id.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.post(path, headers=headers)
        print(response.data)
        self.assertEqual(201, response.status_code)
        response_post_json = json.loads(response.data)
        
        # Retrieves the version_id from the response
        revision_id = response_post_json["id"]

        # Ensures that getting the version has:
        # - PRIVATE visibility (since it's unpublished)
        # - WRITE access_type if the header is from the owner
        # - The response matches the one from the POST request
        path = f"/dp/v1/collections/{revision_id}"
        response = self.app.get(path, headers=headers)
        self.assertEqual(200, response.status_code)
        response_json = json.loads(response.data)
        self.assertEqual("PRIVATE", response_json["visibility"])
        self.assertEqual("WRITE", response_json["access_type"])
        self.assertEqual(response_post_json, response_json)

        # If no auth is passed, the collection should be returned with access_type=READ
        path = f"/dp/v1/collections/{revision_id}"
        response = self.app.get(path)
        self.assertEqual(200, response.status_code)
        response_json = json.loads(response.data)
        self.assertEqual("READ", response_json["access_type"])

    # ✅
    def test__start_revision_of_a_collection_w_links__201(self):
        """
        Start a revision of a collection with links.
        """

        # Create a published collection with 2 datasets and 2 links
        links = [
            Link("Link 1", "OTHER", "http://link.good"),
            Link("DOI Link", "DOI", "http://doi.org/10.1016"),
        ]
        published_collection = self.generate_published_collection(links=links, add_datasets=2)

        # Starts a revision
        path = f"/dp/v1/collections/{published_collection.collection_id.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.post(path, headers=headers)
        self.assertEqual(201, response.status_code)
        response_post_json = json.loads(response.data)

        # Verify that the links are in the revision
        self.assertIsNotNone(response_post_json["links"])
        self.assertCountEqual(
            [link["link_name"] for link in response_post_json["links"]],
            ["Link 1", "DOI Link"] 
        )

    # ✅
    def test__revision__403(self):
        """
        Starting a revision on a revision.
        """
        published_collection = self.generate_published_collection(add_datasets=2)
        test_url = f"/dp/v1/collections/{published_collection.collection_id.id}"

        # Start a revision
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)

        # Try to start a revision again
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)

    # ✅
    def test__revision_nonexistent__403(self):
        """
        Start a revision on a non-existing collection.
        """
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.post("/dp/v1/collections/random", headers=headers)
        self.assertEqual(403, response.status_code)

    # ✅
    def test__revision_not_owner__403(self):
        """
        Start a revision on a collection as a non-owner.
        """
        published_collection = self.generate_published_collection(owner="someone_else", add_datasets=2)
        test_url = f"/dp/v1/collections/{published_collection.collection_id.id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)

class TestDeleteRevision(NewBaseTest):
    """Test case for deleting a collection or datasets under revision."""

    def _create_revision(self, collection_id: str, user="owner"):
        """
        Convenience method for creating a revision on the fly
        """
        path = f"/dp/v1/collections/{collection_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token(user)}
        response = self.app.post(path, headers=headers)
        response_post_json = json.loads(response.data)
        print(response_post_json)
        revision_id = response_post_json["id"]
        return revision_id

    # ✅
    def test__revision_deleted__204(self):
        """
        Delete a collection under revision.
        """
        published_collection = self.generate_published_collection(add_datasets=2)
        revision_id = self._create_revision(published_collection.collection_id.id)

        # Delete the revision
        path = f"/dp/v1/collections/{revision_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        resp = self.app.delete(path, headers=headers)
        self.assertEqual(204, resp.status_code)

        # Cannot get the revision
        path = f"/dp/v1/collections/{revision_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        resp = self.app.get(path, headers=headers)
        self.assertEqual(403, resp.status_code)

        # Verify that the published collection still has two datasets
        path = f"/dp/v1/collections/{published_collection.collection_id.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        resp = self.app.get(path, headers=headers)
        self.assertEqual(200, resp.status_code)
        resp_json = json.loads(resp.data)
        datasets = resp_json["datasets"]
        self.assertEqual(2, len(datasets))

    # ✅
    def test__delete_revision_unauth__403(self):
        """
        An unauthorized user should not be able to delete a collection revision
        """
        published_collection = self.generate_published_collection(owner ="someone_else", add_datasets=2)
        revision_id = self._create_revision(published_collection.collection_id.id, user="super")

        # Delete the revision
        path = f"/dp/v1/collections/{revision_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        resp = self.app.delete(path, headers=headers)
        self.assertEqual(403, resp.status_code)

        
# Those were the previous revision tests - mostly covered by the business logic layer now
# A few cases are good, but we don't need to test every single case here
class TestPublishRevision(NewBaseTest):
    """Test case for publishing a revision."""

    def setUp(self):
        super().setUp()
        self.base_path = "/dp/v1/collections"
        self.mock_timestamp = datetime(2000, 12, 25, 0, 0)
        # TODO: header pattern
        self.headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

    # 💛 passes, but assertions need to be improved
    def test__publish_revision_with_new_dataset__OK(self):
        """
        Publish a revision with new datasets.
        """

        unpublished_collection = self.generate_unpublished_collection(add_datasets=1)

        path = f"{self.base_path}/{unpublished_collection.version_id.id}/publish"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        body = {"data_submission_policy_version": "1.0"} # TODO: still in use?
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        
        self.assertEqual(202, response.status_code)
        self.assertDictEqual(
            {"collection_id": unpublished_collection.collection_id.id, "visibility": "PUBLIC"}, 
            json.loads(response.data)
        )

        # Check GET collection/<collection_id>
        path = f"{self.base_path}/{unpublished_collection.collection_id.id}"
        response = self.app.get(path, headers=headers)
        self.assertEqual(200, response.status_code)

        response_json = json.loads(response.data)
        self.assertEqual("PUBLIC", response_json["visibility"])
        self.assertEqual(unpublished_collection.collection_id.id, response_json["id"])

        # Check published_at and revised_at
        # Collection: Only revised_at should be updated
        self.assertIsNotNone(response_json.get("published_at"))
        # TODO: revised_at
        # TODO: timestamp mocking
        # self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(response_json["revised_at"]))

        # Datasets: Only the newly added dataset should have published_at updated
        # TODO: dataset published_at
        # for dataset in response_json["datasets"]:
        #     if dataset["id"] == new_dataset_id:
        #         self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(dataset["published_at"]))
        #     else:
        #         self.assertIsNone(dataset.get("published_at"))
        #     self.assertIsNone(dataset.get("revised_at"))

    # 💛 passes, but assertions need to be improved
    def test__publish_revision_with_removed_datasets__OK(self):
        """
        Publish a revision with delete datasets.
        """
        unpublished_collection = self.generate_unpublished_collection(add_datasets=2)
        dataset_to_delete = unpublished_collection.datasets[0]

        # Delete a dataset under revision
        self.app.delete(f"/dp/v1/datasets/{dataset_to_delete.version_id.id}", headers=self.headers)

        # Publish the revision with the deleted dataset
        path = f"{self.base_path}/{unpublished_collection.version_id.id}/publish"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        body = {"data_submission_policy_version": "1.0"} # TODO: still in use?
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)

        response_json = json.loads(response.data)
        print(response.data)
        id = response_json["collection_id"]

        # Checks that the published revision exists and only has one dataset
        published_collection = self.business_logic.get_published_collection_version(CollectionId(id))
        self.assertIsNotNone(published_collection)
        if published_collection:
            self.assertIsNotNone(published_collection.published_at)
            self.assertEqual(published_collection.collection_id, unpublished_collection.collection_id)
            self.assertEqual(1, len(published_collection.datasets))

        # TODO: published_at, revised_at

    # ✅
    def test__publish_revision_with_all_removed_datasets__409(self):
        """
        Unable to publish a revision with no datasets.
        """
        unpublished_collection = self.generate_unpublished_collection(add_datasets=1)
        dataset_to_delete = unpublished_collection.datasets[0]

        # Delete a dataset under revision
        self.app.delete(f"/dp/v1/datasets/{dataset_to_delete.version_id.id}", headers=self.headers)

        # Publish the revision with the deleted dataset
        path = f"{self.base_path}/{unpublished_collection.version_id.id}/publish"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        body = {"data_submission_policy_version": "1.0"} # TODO: still in use?
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(409, response.status_code)

    # 💛 Returns 401 instead of 403, verify current behavior
    def test__publish_revision_as_non_autorized__403(self):
        """
        Publish a collection as a non authenticated user returns 403
        """
        collection = self.generate_unpublished_collection(add_datasets=1)
        path = f"{self.base_path}/{collection.version_id.id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    # ✅
    def test__publish_revision_as_not_owner__403(self):
        """
        Publish a collection as a non-owner returns 403
        """
        collection = self.generate_unpublished_collection(add_datasets=1, owner="someone_else")
        path = f"{self.base_path}/{collection.version_id.id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        response = self.app.post(path, headers=self.headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    # ✅
    def test__publish_revision_bad_id__403(self):
        """
        Publish a collection with a bad uuid (non existant) returns 403
        """
        collection_id = "bad_id"
        body = {"data_submission_policy_version": "1.0"}
        path = f"{self.base_path}/{collection_id}/publish"
        response = self.app.post(path, headers=self.headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    # ✅
    def test__can_publish_owned_collection(self):
        collection = self.generate_unpublished_collection(add_datasets=1)
        path = f"{self.base_path}/{collection.version_id.id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("owner")}
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)

    # ✅
    def test__can_publish_non_owned_collection_as_super_curator(self):
        collection = self.generate_unpublished_collection(add_datasets=1, owner="someone else")
        path = f"{self.base_path}/{collection.version_id.id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("super")}
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)

    # The following tests are good to have, but they're essentially business logic tests. We can purge them for now

    # def test__publish_revision_with_updated_datasets__OK(self):
    #     """ "Publish a revision with refreshed datasets."""
    #     self.refresh_datasets()

    #     response_json = self.publish_collection(self.rev_collection)
    #     self.verify_datasets(response_json, {ds.id for ds in self.pub_collection.datasets})

    #     # Check published_at and revised_at
    #     # Collection: Only revised_at should be updated
    #     self.assertIsNone(response_json.get("published_at"))
    #     self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(response_json["revised_at"]))

    #     # Datatsets: None should be updated
    #     for dataset in response_json["datasets"]:
    #         self.assertIsNone(dataset.get("published_at"))
    #         self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(dataset["revised_at"]))

    # def test__publish_revision_with_collection_info_updated__201(self):
    #     """Publish a revision with collection detail changes."""
    #     expected_body = self.update_revision_details()
    #     pub_s3_objects, _ = self.get_s3_objects_from_collections()

    #     # Published revision with collection details updated
    #     response_json = self.publish_collection(self.rev_collection)
    #     self.assertPublishedCollectionOK(expected_body, pub_s3_objects)

    #     self.verify_datasets(response_json, {ds.id for ds in self.pub_collection.datasets})

    #     # Check published_at and revised_at
    #     # Collection: None should be updated
    #     self.assertIsNone(response_json.get("published_at"))
    #     self.assertIsNone(response_json.get("revised_at"))

    #     # Datasets: None should be updated
    #     for dataset in response_json["datasets"]:
    #         self.assertIsNone(dataset.get("published_at"))
    #         self.assertIsNone(dataset.get("revised_at"))

    # def test__publish_revision_with_collection_info_updated_and_refreshed_datasets__201(self):
    #     """Publish a revision with collection detail changes and refreshed datasets."""
    #     expected_body = self.update_revision_details()
    #     self.refresh_datasets()
    #     _, rev_s3_objects = self.get_s3_objects_from_collections()
    #     # Published revision with collection details updated, new dataset, and refreshed datasets
    #     response_json = self.publish_collection(self.rev_collection)
    #     self.verify_datasets(response_json, {ds.id for ds in self.pub_collection.datasets})
    #     self.assertPublishedCollectionOK(expected_body, rev_s3_objects)

    #     # Check published_at and revised_at
    #     # Collection: Only revised_at should be updated
    #     self.assertIsNone(response_json.get("published_at"))
    #     self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(response_json["revised_at"]))

    #     # Datatsets: None should be updated
    #     for dataset in response_json["datasets"]:
    #         self.assertIsNone(dataset.get("published_at"))
    #         self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(dataset["revised_at"]))

    # def test__publish_revision_with_collection_info_updated_and_new_datasets__201(self):
    #     """Publish a revision with collection detail changes and new datasets."""
    #     expected_body = self.update_revision_details()

    #     # add new dataset
    #     new_dataset_id = self.generate_dataset_with_s3_resources(self.session, collection_id=self.rev_collection.id).id
    #     dataset_ids = {ds.id for ds in self.pub_collection.datasets}
    #     dataset_ids.add(new_dataset_id)

    #     # get revision datasets
    #     _, rev_s3_objects = self.get_s3_objects_from_collections()

    #     # Published revision with collection details updated, new dataset, and refreshed datasets
    #     response_json = self.publish_collection(self.rev_collection)
    #     self.assertPublishedCollectionOK(expected_body, rev_s3_objects)

    #     # Check published_at and revised_at
    #     # Collection: Only revised_at should be updated
    #     self.assertIsNone(response_json.get("published_at"))
    #     self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(response_json["revised_at"]))

    #     # Datasets: Only the newly added dataset should have published_at updated
    #     for dataset in response_json["datasets"]:
    #         if dataset["id"] == new_dataset_id:
    #             self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(dataset["published_at"]))
    #         else:
    #             self.assertIsNone(dataset.get("published_at"))
    #         self.assertIsNone(dataset.get("revised_at"))

    # def test__publish_revision_with_collection_info_updated_new_and_refreshed_datasets__201(self):
    #     """Publish a revision with collection detail changes, new datasets, and refreshed datasets."""
    #     expected_body = self.update_revision_details()
    #     self.refresh_datasets()

    #     # add new dataset
    #     new_dataset_id = self.generate_dataset_with_s3_resources(self.session, collection_id=self.rev_collection.id).id
    #     dataset_ids = {ds.id for ds in self.pub_collection.datasets}
    #     dataset_ids.add(new_dataset_id)

    #     # get revision datasets
    #     _, rev_s3_objects = self.get_s3_objects_from_collections()

    #     # Published revision with collection details updated, new dataset, and refreshed datasets
    #     response_json = self.publish_collection(self.rev_collection)
    #     self.assertPublishedCollectionOK(expected_body, rev_s3_objects)

    #     # Check published_at and revised_at
    #     # Collection: Only revised_at should be updated
    #     self.assertIsNone(response_json.get("published_at"))
    #     self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(response_json["revised_at"]))

    #     # Datasets: Only the newly added dataset should have published_at updated
    #     for dataset in response_json["datasets"]:
    #         if dataset["id"] == new_dataset_id:
    #             self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(dataset["published_at"]))
    #             self.assertIsNone(dataset.get("revised_at"))
    #         else:
    #             self.assertIsNone(dataset.get("published_at"))
    #             self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(dataset["revised_at"]))

    # def test__publish_revision_and_existing_datasets(self):
    #     """Publish a revision with the same, existing datasets."""
    #     response_json = self.publish_collection(self.rev_collection)
    #     self.verify_datasets(response_json, {ds.id for ds in self.pub_collection.datasets})

    #     # Check published_at and revised_at
    #     # Collection: None should be updated
    #     self.assertIsNone(response_json.get("published_at"))
    #     self.assertIsNone(response_json.get("revised_at"))

    #     # Datasets: None should be updated
    #     for dataset in response_json["datasets"]:
    #         self.assertIsNone(dataset.get("published_at"))
    #         self.assertIsNone(dataset.get("revised_at"))


###### UPLOAD TESTS START HERE ######

class TestCollectionPostUploadLink(NewBaseTest):

    def setUp(self):
        super().setUp()
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        self.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "file.h5ad"))

    # ✅
    def test__link__202(self):
        collection = self.generate_unpublished_collection()
        path = f"/dp/v1/collections/{collection.version_id}/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        body = {"url": self.good_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)
        self.assertIn("dataset_id", json.loads(response.data).keys())

    # ✅
    def test__link_no_auth__401(self):
        collection = self.generate_unpublished_collection()
        path = f"/dp/v1/collections/{collection.version_id}/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(401, response.status_code)

    # ✅
    def test__link_not_owner__403(self):
        collection = self.generate_unpublished_collection(owner="someone else")
        path = f"/dp/v1/collections/{collection.version_id}/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    # 💛 Needs attention
    def test__bad_link__400(self):
        collection = self.generate_unpublished_collection()
        path = f"/dp/v1/collections/{collection.version_id}/upload-links"

        # ✅
        with self.subTest("Unsupported Provider"):
            # Mocks the URI validator so that it return invalid
            self.uri_provider.validate = Mock(return_value=False)
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
            body = {"url": "https://test_url.com"}
            test_url = furl(path=path)
            response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
            self.assertEqual(400, response.status_code)
            self.assertEqual("The dropbox shared link is invalid.", json.loads(response.data)["detail"])

        # 🔴 TODO: requires a finer grained exception management
        with self.subTest("Bad Dropbox link"):
            self.uri_provider.validate = Mock(return_value=True)

            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
            body = {"url": self.dummy_link}
            test_url = furl(path=path)
            response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
            self.assertEqual(400, response.status_code)
            self.assertTrue(
                "'content-disposition' not present in the header." in json.loads(response.data)["detail"]
                or "The URL provided causes an error with Dropbox." == json.loads(response.data)["detail"]
            )

    # ✅
    def test__oversized__413(self):
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(40 * 2**30, "file.h5ad"))
        collection = self.generate_unpublished_collection()
        path = f"/dp/v1/collections/{collection.version_id}/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(413, response.status_code)

    # ✅
    def test__link_fake_collection__403(self):
        path = "/dp/v1/collections/fake/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        body = {"url": self.good_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    # ✅
    def test_link_public_collection__403(self):
        collection = self.generate_published_collection()
        path = f"/dp/v1/collections/{collection.version_id}/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        body = {"url": self.good_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)


class TestCollectionPutUploadLink(NewBaseTest):
    def setUp(self):
        super().setUp()
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        # TODO: headers do not follow the same pattern as the rest of the test. Maybe change?
        self.headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "file.h5ad"))

    def assert_datasets_are_updated(self, dataset1, dataset2):
        """
        Asserts that dataset2 is an updated revision of dataset1
        """
        self.assertIsNotNone(dataset2)
        if dataset2 is not None: #pylance
            # The new dataset should have the same canonical id
            self.assertNotEqual(dataset1.dataset_version_id, dataset2.version_id)
            self.assertEqual(dataset2.dataset_id.id, dataset1.dataset_id)
            # The artifacts for the new datasets should have a different id

    # ✅
    def test__reupload_published_dataset_during_revision__202(self):
        """
        Reupload a published dataset during a revision
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("processing_status", DatasetProcessingStatus.SUCCESS)],
            publish=True
        )

        new_version = self.business_logic.create_collection_version(CollectionId(dataset.collection_id))

        path = f"/dp/v1/collections/{new_version.version_id}/upload-links"
        body = {"url": self.good_link, "id": dataset.dataset_version_id}
        response = self.app.put(path, headers=self.headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)

        new_dataset_id = json.loads(response.data)["dataset_id"]
        new_dataset = self.business_logic.get_dataset_version(DatasetVersionId(new_dataset_id))
        self.assert_datasets_are_updated(dataset, new_dataset)

    # ✅
    @patch("backend.corpora.common.upload.start_upload_sfn")
    def test__reupload_unpublished_dataset__202(self):
        """
        Reuploads an unpublished dataset
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("processing_status", DatasetProcessingStatus.SUCCESS)],
        )

        path = f"/dp/v1/collections/{dataset.collection_version_id}/upload-links"
        body = {"url": self.good_link, "id": dataset.dataset_version_id}

        response = self.app.put(path, headers=self.headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)
        
        new_dataset_id = json.loads(response.data)["dataset_id"]
        new_dataset = self.business_logic.get_dataset_version(DatasetVersionId(new_dataset_id))
        self.assertNotEqual(new_dataset_id, dataset.dataset_version_id)
        self.assert_datasets_are_updated(dataset, new_dataset)

    # ✅
    def test__reupload_public_dataset__403(self):
        """cannot reupload a public published dataset"""
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("processing_status", DatasetProcessingStatus.SUCCESS)],
            publish=True,
        )
        path = f"/dp/v1/collections/{dataset.collection_version_id}/upload-links"
        body = {"url": self.good_link, "id": dataset.dataset_version_id}

        response = self.app.put(path, headers=self.headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)
        # TODO: possibly, assert that no new version has been created, but not too useful

    # ✅
    def test__reupload_while_processing_dataset__405(self):
        """cannot reupload a dataset that is pending"""
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate("processing_status", DatasetProcessingStatus.PENDING)],
        )
        path = f"/dp/v1/collections/{dataset.collection_version_id}/upload-links"
        body = {"url": self.good_link, "id": dataset.dataset_version_id}
        
        response = self.app.put(path, headers=self.headers, data=json.dumps(body))
        print(response.data)
        self.assertEqual(405, response.status_code)

    # ✅
    def test__reupload_dataset_not_owner__403(self):
        dataset = self.generate_dataset(
            owner="someone else",
            statuses=[DatasetStatusUpdate("processing_status", DatasetProcessingStatus.SUCCESS)],
        )
        
        path = f"/dp/v1/collections/{dataset.collection_version_id}/upload-links"
        body = {"url": self.good_link, "id": dataset.dataset_version_id}

        response = self.app.put(path, headers=self.headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    # ✅
    def test__dataset_not_in_collection__404(self):
        """
        Trying to publish a dataset belonging to a different collection raises a 404
        """
        dataset = self.generate_dataset(
            owner="someone else",
            statuses=[DatasetStatusUpdate("processing_status", DatasetProcessingStatus.SUCCESS)],
        )
        collection = self.generate_unpublished_collection()
        
        path = f"/dp/v1/collections/{collection.version_id.id}/upload-links"
        body = {"url": self.good_link, "id": dataset.dataset_version_id}

        response = self.app.put(path, headers=self.headers, data=json.dumps(body))
        self.assertEqual(404, response.status_code)


class TestCollectionUploadLinkCurators(NewBaseTest):
    def setUp(self):
        super().setUp()
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "file.h5ad"))

    # ✅
    def test__can_upload_dataset_to_non_owned_collection_as_super_curator(self):
        """
        A super curator can upload a dataset to a non-owned collection
        """
        collection = self.generate_unpublished_collection(
            owner="someone else"
        )
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("super")}
        path = f"/dp/v1/collections/{collection.version_id}/upload-links"
        body = {"url": self.good_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        print(response.data)
        self.assertEqual(202, response.status_code)
        self.assertIn("dataset_id", json.loads(response.data).keys())

    # ✅
    def test__can_reupload_dataset_as_super_curator(self):
        """
        A super curator can reupload a dataset to a non-owned collection
        """
        dataset = self.generate_dataset(
            owner="someone else",
            statuses=[DatasetStatusUpdate("processing_status", DatasetProcessingStatus.SUCCESS)],
        )

        path = f"/dp/v1/collections/{dataset.collection_version_id}/upload-links"
        body = {"url": self.good_link, "id": dataset.dataset_version_id}

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("super")}
        response = self.app.put(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)
