import itertools
import json
import unittest
from datetime import datetime
from unittest.mock import patch

from furl import furl

from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
    UploadStatus,
    generate_id,
    ProjectLinkType,
)
from backend.corpora.common.entities import Collection
from backend.corpora.common.providers.crossref_provider import CrossrefDOINotFoundException, CrossrefFetchException
from backend.corpora.lambdas.api.v1.collection import verify_collection_body
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest, get_cxguser_token
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


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


class TestCollection(BaseAuthAPITest):
    def validate_collections_response_structure(self, body):
        self.assertIn("collections", body)
        self.assertTrue(all(k in ["collections", "from_date", "to_date"] for k in body))

        for collection in body["collections"]:
            self.assertListEqual(sorted(collection.keys()), ["created_at", "id", "visibility"])
            self.assertEqual(collection["visibility"], "PUBLIC")
            self.assertGreaterEqual(datetime.fromtimestamp(collection["created_at"]).year, 1969)

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

    def test__get_collection_id__ok(self):
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
                            "user_submitted": True,
                        },
                        {
                            "dataset_id": "test_dataset_id",
                            "filename": "test_filename",
                            "filetype": "CXG",
                            "id": "test_dataset_artifact_id_cxg",
                            "s3_uri": "s3://bogus-bucket/test_s3_uri.h5ad",
                            "user_submitted": True,
                        },
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
                    "organism": [{"label": "test_organism", "ontology_term_id": "test_obo"}],
                    "collection_id": "test_collection_id",
                    "cell_type": [{"label": "test_cell_type", "ontology_term_id": "test_opo"}],
                    "x_normalization": "test_x_normalization",
                    "x_approximate_distribution": "NORMAL",
                    "batch_condition": ["batchA", "batchB"],
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
                        "h5ad_status": "NA",
                        "rds_status": "NA",
                        "cxg_status": "NA",
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
                    "linked_datasets": [],
                    "description": "this is a geneset",
                    "id": "test_geneset",
                    "name": "test_geneset",
                },
                {
                    "collection_id": "test_collection_id",
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
            "contact_name": "Some Body",
            "curator_name": "",
            "contact_email": "somebody@chanzuckerberg.com",
            "data_submission_policy_version": "0",
        }

        with self.subTest("auth cookie"):
            expected_body["access_type"] = "WRITE"
            test_url = furl(path="/dp/v1/collections/test_collection_id", query_params=dict(visibility="PUBLIC"))
            cxguser_cookie = get_cxguser_token()
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
                "curator_name": "",
                "data_submission_policy_version": "0",
                "datasets": [
                    {
                        "collection_id": collection.id,
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
                    headers["Cookie"] = get_cxguser_token()
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

    def test__get_collection_id__403_not_found(self):
        """Verify the test collection exists and the expected fields exist."""
        test_url = furl(path="/dp/v1/collections/AAAA-BBBB-CCCC-DDDD", query_params=dict(visibility="PUBLIC"))
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, response.status_code)

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
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)

        # Add curator_name
        data["curator_name"] = "john smith"
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)

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
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        collection = Collection.get_collection(
            self.session, collection_id, CollectionVisibility.PRIVATE.name, include_tombstones=True
        )
        self.assertEquals(collection.get_doi(), "https://doi.org/10.1016/foo")

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
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(400, response.status_code)

    @patch("backend.corpora.common.providers.crossref_provider.CrossrefProvider.fetch_metadata")
    def test__post_collection_rejects_doi_not_in_crossref(self, mock_provider):
        mock_provider.side_effect = CrossrefDOINotFoundException("Mocked CrossrefDOINotFoundException")
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
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(400, response.status_code)
        error_payload = json.loads(response.data)
        self.assertEqual(error_payload["detail"][0], {"link_type": "DOI", "reason": "DOI cannot be found on Crossref"})

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
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(400, response.status_code)
        error_payload = json.loads(response.data)
        self.assertEqual([{"link_type": "DOI", "reason": "Invalid DOI"}], error_payload["detail"])

    @patch("backend.corpora.common.providers.crossref_provider.CrossrefProvider.fetch_metadata")
    def test__post_collection_ignores_metadata_if_crossref_exception(self, mock_provider):
        mock_provider.side_effect = CrossrefFetchException("Mocked CrossrefFetchException")
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
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        collection = Collection.get_collection(
            self.session, collection_id, CollectionVisibility.PRIVATE.name, include_tombstones=True
        )
        self.assertIsNone(collection.publisher_metadata)

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
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        collection = Collection.get_collection(
            self.session, collection_id, CollectionVisibility.PRIVATE.name, include_tombstones=True
        )
        self.assertIsNone(collection.publisher_metadata)

    @patch("backend.corpora.common.providers.crossref_provider.CrossrefProvider.fetch_metadata")
    def test__post_collection_adds_publisher_metadata(self, mock_provider):

        mock_provider.return_value = generate_mock_publisher_metadata()

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
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        collection = Collection.get_collection(
            self.session, collection_id, CollectionVisibility.PRIVATE.name, include_tombstones=True
        )
        self.assertEqual(collection.publisher_metadata, generate_mock_publisher_metadata())

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
                        "Cookie": get_cxguser_token(),
                    },
                    data=json.dumps(data),
                )
                self.assertEqual(400, response.status_code)

    def test__post_collection_fails_if_data_missing(self):
        test_url = furl(path="/dp/v1/collections")
        data = json.dumps({"name": "bkjbjbjmbjm"})
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()},
            data=data,
        )
        self.assertEqual(400, response.status_code)

    def test__can_retrieve_created_collection(self):
        test_url = furl(path="/dp/v1/collections")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
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
        test_url.add(query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers=no_cookie_headers)
        self.assertEqual("READ", json.loads(response.data)["access_type"])

    def test__create_collection_strip_string_fields(self):
        test_url = furl(path="/dp/v1/collections")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
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
        revision_not_owned = self.generate_collection(
            self.session,
            visibility=CollectionVisibility.PRIVATE.name,
            owner="someone else",
            revision_of=public_not_owned,
        ).id
        revision_owned = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id", revision_of=public_owned
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
            self.assertNotIn(revision_owned, ids)

        with self.subTest("auth"):
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
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
            self.assertNotIn(revision_not_owned, ids)
            self.assertIn(revision_owned, ids)
            self.assertTrue(
                [collection for collection in collections if collection.get("revision_of") == public_owned][0]
            )

    @unittest.skip("Flaky")  # TODO @mdunitz ticket 2223
    def test__get_all_collections_for_index(self):
        collection = self.generate_collection(
            self.session,
            visibility=CollectionVisibility.PUBLIC.name,
            owner="test_user_id",
            publisher_metadata=generate_mock_publisher_metadata(),
            published_at=datetime.now(),
            revised_at=datetime.now(),
        )

        collection_id_private = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        ).id

        collection_id_tombstone = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, tombstone=True, owner="test_user_id"
        ).id

        test_url = furl(path="/dp/v1/collections/index")
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [collection["id"] for collection in body]
        self.assertIn(collection.id, ids)
        self.assertNotIn(collection_id_private, ids)
        self.assertNotIn(collection_id_tombstone, ids)

        actual_collection = body[-1]  # last added collection
        self.assertEqual(actual_collection["id"], collection.id)
        self.assertEqual(actual_collection["name"], collection.name)
        self.assertNotIn("description", actual_collection)
        self.assertEqual(actual_collection["published_at"], collection.published_at.timestamp())
        self.assertEqual(actual_collection["revised_at"], collection.revised_at.timestamp())
        self.assertEqual(actual_collection["publisher_metadata"], collection.publisher_metadata)

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
                headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
                response = self.app.post("/dp/v1/collections", headers=headers, data=json.dumps(body))
                self.assertEqual(400, response.status_code)
                for error in expected_errors:
                    self.assertIn(error, response.json["detail"])


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


class TestUpdateCollection(BaseAuthAPITest):
    def test__update_collection__OK(self):
        collection = self.generate_collection(self.session)
        test_fields = [
            "name",
            "description",
            "contact_name",
            "contact_email",
            "links",
        ]
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}

        # Update the collection
        expected_body = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
        }
        data = json.dumps(expected_body)
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        for field in test_fields:
            self.assertEqual(expected_body[field], actual_body[field])

    def test__update_collection__403(self):
        collection = self.generate_collection(self.session, owner="someone else")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        data = json.dumps({"name": "new name"})
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)
        self.assertEqual(403, response.status_code)

    def test__update_collection_links__OK(self):
        collection = self.generate_collection(self.session)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}

        # add links
        links = [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}]
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
            {"link_name": "Link 1", "link_url": "http://link.com", "link_type": "OTHER"},
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

    @patch("backend.corpora.common.providers.crossref_provider.CrossrefProvider.fetch_metadata")
    def test__update_collection_new_doi_updates_metadata(self, mock_provider):
        # The Crossref provider will always return "New Journal"
        mock_provider.return_value = generate_mock_publisher_metadata("New Journal")
        collection = self.generate_collection(
            self.session,
            links=[{"link_name": "Link 1", "link_url": "http://doi.org/123", "link_type": "DOI"}],
            publisher_metadata=generate_mock_publisher_metadata("Old Journal"),
        )
        self.assertEqual("Old Journal", collection.publisher_metadata["journal"])

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.put(
            f"/dp/v1/collections/{collection.id}",
            data=json.dumps({"links": [{"link_name": "Link 1", "link_url": "http://doi.org/456", "link_type": "DOI"}]}),
            headers=headers,
        )
        self.assertEqual(200, response.status_code)

        mock_provider.assert_called_once()

        actual_body = json.loads(response.data)
        self.assertIsNotNone(actual_body["publisher_metadata"])
        self.assertIsNotNone(actual_body["publisher_metadata"]["journal"])
        self.assertEqual("New Journal", actual_body["publisher_metadata"]["journal"])

    @patch("backend.corpora.common.providers.crossref_provider.CrossrefProvider.fetch_metadata")
    def test__update_collection_remove_doi_deletes_metadata(self, mock_provider):
        mock_provider.return_value = generate_mock_publisher_metadata("New Journal")
        collection = self.generate_collection(
            self.session,
            links=[{"link_name": "Link 1", "link_url": "http://doi.org/123", "link_type": "DOI"}],
            publisher_metadata=generate_mock_publisher_metadata("Old Journal"),
        )
        self.assertEqual("Old Journal", collection.publisher_metadata["journal"])

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        # We're passing an empty links object, therefore the DOI is deleted
        response = self.app.put(
            f"/dp/v1/collections/{collection.id}",
            data=json.dumps({"links": []}),
            headers=headers,
        )

        self.assertEqual(200, response.status_code)
        # No Crossref calls should happen
        mock_provider.assert_not_called()

        # The `publisher_metadata` node should not longer be in the collection
        actual_body = json.loads(response.data)
        self.assertNotIn("publisher_metadata", actual_body)

    @patch("backend.corpora.common.providers.crossref_provider.CrossrefProvider.fetch_metadata")
    def test__update_collection_same_doi_does_not_update_metadata(self, mock_provider):
        mock_provider.return_value = generate_mock_publisher_metadata("New Journal")
        collection = self.generate_collection(
            self.session,
            links=[{"link_name": "Link 1", "link_url": "http://doi.org/123", "link_type": "DOI"}],
            publisher_metadata=generate_mock_publisher_metadata("Old Journal"),
        )
        self.assertEqual("Old Journal", collection.publisher_metadata["journal"])

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        # Note that the DOI is the same as the original
        response = self.app.put(
            f"/dp/v1/collections/{collection.id}",
            data=json.dumps({"links": [{"link_name": "Link 1", "link_url": "http://doi.org/123", "link_type": "DOI"}]}),
            headers=headers,
        )

        self.assertEqual(200, response.status_code)
        # No Crossref calls should happen, since the DOI is unmodified
        mock_provider.assert_not_called()

        # The `publisher_metadata` node should exist and be the same
        actual_body = json.loads(response.data)
        self.assertIsNotNone(actual_body["publisher_metadata"])
        self.assertIsNotNone(actual_body["publisher_metadata"]["journal"])
        self.assertEqual("Old Journal", actual_body["publisher_metadata"]["journal"])


class TestCollectionsCurators:
    def test_view_non_owned_private_collection__ok(self):
        # Generate test collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="another_test_user_id"
        )

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        test_url = furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers=headers)

        # This will pass even for non curators.
        # Why are users allowed to view private collections that they don't own?
        self.assertEqual(response.status_code, 200)

        body = json.loads(response.data)
        self.assertEqual(body["access_type"], "WRITE")

    def test_update_non_owned_private_collection_as_super_curator__ok(self):
        # Generate test collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="another_test_user_id"
        )

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token("super")}

        modified_collection = {
            "name": "new name",
            "description": collection.description,
            "contact_name": collection.contact_name,
            "contact_email": collection.contact_email,
            "links": collection.links,
        }
        data = json.dumps(modified_collection)
        response = self.app.put(f"/dp/v1/collections/{collection.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)

    def test_delete_non_owned_private_collection_as_super_curator__ok(self):
        # Generate test collection
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="another_test_user_id"
        )

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token("super")}

        test_url = furl(path=f"/dp/v1/collections/{collection.id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(204, response.status_code)


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
