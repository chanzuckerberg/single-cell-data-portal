import dataclasses
import itertools
import json
from datetime import datetime
from unittest import mock
from unittest.mock import Mock, patch

from furl import furl

from backend.layers.business.entities import DatasetArtifactDownloadData
from backend.layers.common.entities import (
    CollectionId,
    CollectionVersionId,
    DatasetArtifactType,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetVersionId,
    Link,
    OntologyTermId,
)
from backend.layers.thirdparty.crossref_provider import CrossrefDOINotFoundException, CrossrefFetchException
from backend.layers.thirdparty.uri_provider import FileInfo, FileInfoException
from tests.unit.backend.layers.api.fixture import generate_mock_publisher_metadata
from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest
from tests.unit.backend.layers.common.base_test import DatasetArtifactUpdate, DatasetStatusUpdate


class TestCollection(BaseAPIPortalTest):
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
            self.assertIn(expected_id, [p["id"] for p in actual_body["collections"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(None, actual_body.get("from_date"))

    # ✅
    def test__get_collection_id__ok(self):
        """Verify the test collection exists and the expected fields exist."""

        collection = self.generate_published_collection(
            add_datasets=2,
        )

        self.maxDiff = None

        expected_body = {
            "access_type": "WRITE",
            "consortia": ["Consortia 1", "Consortia 2"],
            "contact_email": "john.doe@email.com",
            "contact_name": "john doe",
            "created_at": mock.ANY,
            "curator_name": "Jane Smith",
            "data_submission_policy_version": "1.0",
            "datasets": [
                {
                    "assay": [{"label": "test_assay_label", "ontology_term_id": "test_assay_term_id"}],
                    "batch_condition": ["test_batch_1", "test_batch_2"],
                    "cell_count": 10,
                    "primary_cell_count": 5,
                    "cell_type": [{"label": "test_cell_type_label", "ontology_term_id": "test_cell_type_term_id"}],
                    "collection_id": collection.collection_id.id,
                    "created_at": mock.ANY,
                    "dataset_assets": mock.ANY,
                    "dataset_deployments": [{"url": mock.ANY}],
                    "development_stage": [
                        {"label": "test_development_stage_label", "ontology_term_id": "test_development_stage_term_id"}
                    ],
                    "disease": [{"label": "test_disease_label", "ontology_term_id": "test_disease_term_id"}],
                    "donor_id": ["test_donor_1"],
                    "id": mock.ANY,
                    "is_primary_data": "BOTH",
                    "is_valid": True,
                    "mean_genes_per_cell": 0.5,
                    "name": "test_dataset_name",
                    "organism": [{"label": "test_organism_label", "ontology_term_id": "test_organism_term_id"}],
                    "processing_status": {
                        "created_at": 0,
                        "cxg_status": "NA",
                        "dataset_id": mock.ANY,
                        "h5ad_status": "NA",
                        "id": "NA",
                        "processing_status": "INITIALIZED",
                        "rds_status": "NA",
                        "updated_at": 0,
                        "upload_progress": 1,
                        "upload_status": "WAITING",
                        "validation_status": "NA",
                    },
                    "published": True,
                    "published_at": mock.ANY,
                    "revision": 0,  # NA
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
                    "tombstone": False,
                    "updated_at": mock.ANY,
                    "x_approximate_distribution": "normal",
                },
                {
                    "assay": [{"label": "test_assay_label", "ontology_term_id": "test_assay_term_id"}],
                    "batch_condition": ["test_batch_1", "test_batch_2"],
                    "cell_count": 10,
                    "primary_cell_count": 5,
                    "cell_type": [{"label": "test_cell_type_label", "ontology_term_id": "test_cell_type_term_id"}],
                    "collection_id": collection.collection_id.id,
                    "created_at": mock.ANY,
                    "dataset_assets": mock.ANY,
                    "dataset_deployments": [{"url": mock.ANY}],
                    "development_stage": [
                        {"label": "test_development_stage_label", "ontology_term_id": "test_development_stage_term_id"}
                    ],
                    "disease": [{"label": "test_disease_label", "ontology_term_id": "test_disease_term_id"}],
                    "donor_id": ["test_donor_1"],
                    "id": mock.ANY,
                    "is_primary_data": "BOTH",
                    "is_valid": True,
                    "mean_genes_per_cell": 0.5,
                    "name": "test_dataset_name",
                    "organism": [{"label": "test_organism_label", "ontology_term_id": "test_organism_term_id"}],
                    "processing_status": {
                        "created_at": 0,
                        "cxg_status": "NA",
                        "dataset_id": mock.ANY,
                        "h5ad_status": "NA",
                        "id": "NA",
                        "processing_status": "INITIALIZED",
                        "rds_status": "NA",
                        "updated_at": 0,
                        "upload_progress": 1,
                        "upload_status": "WAITING",
                        "validation_status": "NA",
                    },
                    "published": True,
                    "published_at": mock.ANY,
                    "revision": 0,
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
                    "tombstone": False,
                    "updated_at": mock.ANY,
                    "x_approximate_distribution": "normal",
                },
            ],
            "description": "described",
            "id": mock.ANY,
            "links": [],
            "name": "test_collection",
            "published_at": mock.ANY,
            "updated_at": mock.ANY,
            "visibility": "PUBLIC",
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
        test_cases = itertools.product(authenticated, owns, visibility)

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

    def test__get_collection_id_returns_revision_of_published_collection(self):
        version = self.generate_published_collection()
        revision = self.generate_revision(version.collection_id)
        test_url = furl(path=f"/dp/v1/collections/{revision.version_id}")
        headers = dict(host="localhost")
        headers["Cookie"] = self.get_cxguser_token()
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 200)

        body = json.loads(response.data)
        self.assertEqual(body["visibility"], "PRIVATE")
        self.assertEqual(body["access_type"], "WRITE")

        # check revising_in is set if the collection has a revision and the
        # user is logged in and has write access.
        test_url = furl(path=f"/dp/v1/collections/{version.version_id}")
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 200)
        body = json.loads(response.data)
        self.assertEqual(body["visibility"], "PUBLIC")
        self.assertEqual(body["access_type"], "WRITE")
        self.assertEqual(body["revising_in"], revision.version_id.id)

        # check revising_in is not set if the collection has a revision and the
        # user is not logged in.
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(response.status_code, 200)
        body = json.loads(response.data)
        self.assertEqual(body["visibility"], "PUBLIC")
        self.assertEqual(body["access_type"], "READ")
        self.assertNotIn("revising_in", body)

    def test__get_collection_id_retrieves_published_version_by_collection_id(self):
        """
        GET /collections/:collection_id retrieves the published version given a canonical collection_id
        """
        version = self.generate_published_collection()
        test_url = furl(path=f"/dp/v1/collections/{version.collection_id}")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        self.assertEqual(body["visibility"], "PUBLIC")

    def test__get_collection_id_retrieves_published_version_by_collection_id_if_revision(self):
        """
        When there is a revision of a published collection, GET /collections/:collection_id retrieves:
        1. the published version given the canonical collection_id
        2. the revision if given the version_id
        """
        version = self.generate_published_collection()
        revision = self.business_logic.create_collection_version(version.collection_id)

        test_url = furl(path=f"/dp/v1/collections/{revision.collection_id}")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        self.assertEqual(body["visibility"], "PUBLIC")

        test_url = furl(path=f"/dp/v1/collections/{revision.version_id}")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        self.assertEqual(body["visibility"], "PRIVATE")

    def test__get_collection_id_retrieves_unpublished_version(self):
        """
        If a collection is unpublished, GET /collections/:collection_id retrieves the unpublished version
        when passed the canonical collection_id
        """
        version = self.generate_unpublished_collection()

        test_url = furl(path=f"/dp/v1/collections/{version.collection_id}")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        self.assertEqual(body["visibility"], "PRIVATE")

    def test__get_collection_id__403_not_found(self):
        """
        GET /collections/:collection_id returns 403 if an invalid is specified
        """
        fake_id = CollectionId()
        test_url = furl(path=f"/dp/v1/collections/{fake_id}", query_params=dict(visibility="PUBLIC"))
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
            "curator_name": "Curator Name",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
            "consortia": ["Consortia 1"],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)

        # Check that the collection_id is the canonical collection ID
        collection_id = response.json["collection_id"]
        version = self.business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
        self.assertEqual(version.collection_id.id, collection_id)

        # Add curator_name
        data["curator_name"] = "john smith"
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)

    def test__post_collection_rejects_invalid_consortia(self):
        test_url = furl(path="/dp/v1/collections")

        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "curator_name": "Curator Name",
            "links": [],
            "consortia": ["Invalid Consortia"],
        }

        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )

        self.assertEqual(400, response.status_code)

    def test__post_collection_sorts_consortia(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "curator_name": "Curator Name",
            "links": [
                {"link_name": "DOI Link", "link_url": "10.1016/foo", "link_type": "DOI"},
            ],
            "consortia": ["Consortia 3", "Consortia 1"],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        collection = self.business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
        self.assertEqual(collection.metadata.consortia, sorted(data["consortia"]))

    # ✅
    def test__post_collection_normalizes_doi(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "curator_name": "Curator Name",
            "links": [
                {"link_name": "DOI Link", "link_url": "10.1016/foo", "link_type": "DOI"},
            ],
            "consortia": ["Consortia 1"],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        # TODO: this endpoint should also return `version_id`
        collection = self.business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
        doi = next(link.uri for link in collection.metadata.links if link.type == "DOI")  # TODO: careful
        self.assertEquals(doi, "https://doi.org/10.1016/foo")

    # ✅
    def test__post_collection_rejects_two_dois(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "curator_name": "Curator Name",
            "links": [
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
            ],
            "consortia": ["Consortia 1"],
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
            "curator_name": "Curator Name",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
            "consortia": ["Consortia 1"],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(400, response.status_code)
        error_payload = json.loads(response.data)
        self.assertEqual(error_payload["detail"][0], {"link_type": "doi", "reason": "DOI cannot be found on Crossref"})

    # ✅
    def test__post_collection_rejects_invalid_doi(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "curator_name": "Curator Name",
            "links": [{"link_name": "DOI Link", "link_url": "invalid/doi", "link_type": "DOI"}],
            "consortia": ["Consortia 1"],
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
            "curator_name": "Curator Name",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
            "consortia": ["Consortia 1"],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        # TODO: this endpoint should also return `version_id`
        collection = self.business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
        self.assertIsNone(collection.publisher_metadata)

    # ✅
    def test__post_collection_ignores_metadata_if_no_doi(self):
        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "curator_name": "Curator Name",
            "consortia": ["Consortia 1"],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        # TODO: this endpoint should also return `version_id`
        collection = self.business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
        self.assertIsNone(collection.publisher_metadata)

    # ✅
    def test__post_collection_adds_publisher_metadata(self):
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata())

        test_url = furl(path="/dp/v1/collections")
        data = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "curator_name": "Curator Name",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
            "consortia": ["Consortia 1"],
        }
        json_data = json.dumps(data)
        response = self.app.post(
            test_url.url,
            headers={"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()},
            data=json_data,
        )
        self.assertEqual(201, response.status_code)
        collection_id = json.loads(response.data)["collection_id"]
        # TODO: this endpoint should also return `version_id`
        collection = self.business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
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
        self.assertEqual(body["contact_name"], data["contact_name"])
        self.assertEqual(body["contact_email"], data["contact_email"])
        self.assertEqual(body["consortia"], [])

        # test that non owners only have read access
        no_cookie_headers = {"host": "localhost", "Content-Type": "application/json"}
        test_url = furl(path=f"/dp/v1/collections/{collection_id}")
        response = self.app.get(test_url.url, headers=no_cookie_headers)
        self.assertEqual("READ", json.loads(response.data)["access_type"])

        # test that owners have write access
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual("WRITE", json.loads(response.data)["access_type"])

        # test that super curators have write access
        super_curator_headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token("super"),
        }
        response = self.app.get(test_url.url, headers=super_curator_headers)
        self.assertEqual("WRITE", json.loads(response.data)["access_type"])

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
            "consortia": ["Consortia 1   "],
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
        self.assertEqual(body["contact_name"], data["contact_name"].strip())
        self.assertEqual(body["contact_email"], body["contact_email"].strip())
        self.assertEqual(body["data_submission_policy_version"], body["data_submission_policy_version"].strip())
        self.assertEqual(body["consortia"], ["Consortia 1"])

        for link in body["links"]:
            self.assertEqual(link["link_url"], link["link_url"].strip())

    def test__list_collection__check_owner__no_auth(self):
        # Generate test collection
        public_owned = self.generate_published_collection(owner="test_user_id")
        private_owned = self.generate_unpublished_collection(owner="test_user_id")
        public_not_owned = self.generate_published_collection(owner="someone else")
        private_not_owned = self.generate_unpublished_collection(owner="someone else")

        revision_not_owned = self.business_logic.create_collection_version(public_not_owned.collection_id)
        revision_owned = self.business_logic.create_collection_version(public_owned.collection_id)

        path = "/dp/v1/collections"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(path, headers=headers)
        self.assertEqual(200, response.status_code)
        result = json.loads(response.data)
        collections = result.get("collections")
        self.assertIsNotNone(collections)
        ids = [collection.get("id") for collection in collections]

        self.assertIn(public_owned.collection_id.id, ids)
        self.assertNotIn(public_owned.version_id.id, ids)
        self.assertIn(public_not_owned.collection_id.id, ids)
        self.assertNotIn(public_not_owned.version_id.id, ids)
        self.assertNotIn(private_owned.collection_id.id, ids)
        self.assertNotIn(private_owned.version_id.id, ids)
        self.assertNotIn(private_not_owned.collection_id.id, ids)
        self.assertNotIn(private_not_owned.version_id.id, ids)
        self.assertNotIn(revision_owned.version_id.id, ids)
        self.assertNotIn(revision_not_owned.version_id.id, ids)

    def test__list_collection__check_owner__auth(self):
        # Generate test collection
        public_owned = self.generate_published_collection(owner="test_user_id")
        private_owned = self.generate_unpublished_collection(owner="test_user_id")
        public_not_owned = self.generate_published_collection(owner="someone else")
        private_not_owned = self.generate_unpublished_collection(owner="someone else")

        revision_not_owned = self.business_logic.create_collection_version(public_not_owned.collection_id)
        revision_owned = self.business_logic.create_collection_version(public_owned.collection_id)

        path = "/dp/v1/collections"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.get(path, headers=headers)
        self.assertEqual(200, response.status_code)
        result = json.loads(response.data)
        collections = result.get("collections")
        self.assertIsNotNone(collections)
        ids = [collection.get("id") for collection in collections]
        revision_ids = [collection.get("revision_of") for collection in collections if collection.get("revision_of")]
        self.assertIn(public_owned.collection_id.id, ids)
        self.assertNotIn(public_owned.version_id.id, ids)
        self.assertIn(public_not_owned.collection_id.id, ids)
        self.assertNotIn(public_not_owned.version_id.id, ids)
        self.assertIn(private_owned.collection_id.id, ids)
        self.assertNotIn(private_owned.version_id.id, ids)
        self.assertNotIn(private_not_owned.collection_id.id, ids)
        self.assertNotIn(private_not_owned.version_id.id, ids)
        self.assertIn(revision_owned.version_id.id, ids)
        self.assertNotIn(revision_not_owned.version_id.id, ids)
        self.assertIn(public_owned.collection_id.id, revision_ids)

    # ✅
    def test__get_all_collections_for_index(self):
        """
        The `collections/index` endpoint should only return public collections
        """

        collection = self.generate_published_collection()
        collection_to_tombstone = self.generate_published_collection()
        private_collection = self.generate_unpublished_collection()

        self.business_logic.tombstone_collection(collection_to_tombstone.collection_id)

        test_url = furl(path="/dp/v1/collections/index")
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [collection["id"] for collection in body]
        self.assertIn(collection.collection_id.id, ids)
        self.assertNotIn(private_collection.collection_id.id, ids)
        self.assertNotIn(private_collection.version_id.id, ids)
        self.assertNotIn(collection_to_tombstone.collection_id.id, ids)
        self.assertNotIn(collection_to_tombstone.version_id.id, ids)

        actual_collection = body[-1]  # last added collection
        self.assertEqual(actual_collection["id"], collection.collection_id.id)
        self.assertEqual(actual_collection["name"], collection.metadata.name)
        self.assertNotIn("description", actual_collection)
        # Both `published_at` and `revised_at` should point to the same timestamp
        self.assertEqual(actual_collection["published_at"], collection.published_at.timestamp())
        self.assertEqual(actual_collection["revised_at"], collection.published_at.timestamp())

        # test that the owner of a private collection can not see the private collection
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [collection["id"] for collection in body]
        self.assertIn(collection.collection_id.id, ids)
        self.assertNotIn(private_collection.collection_id.id, ids)
        self.assertNotIn(private_collection.version_id.id, ids)
        self.assertNotIn(collection_to_tombstone.collection_id.id, ids)
        self.assertNotIn(collection_to_tombstone.version_id.id, ids)

        # test that super curators can not see the private collection
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("super")}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [collection["id"] for collection in body]
        self.assertIn(collection.collection_id.id, ids)
        self.assertNotIn(private_collection.collection_id.id, ids)
        self.assertNotIn(private_collection.version_id.id, ids)
        self.assertNotIn(collection_to_tombstone.collection_id.id, ids)
        self.assertNotIn(collection_to_tombstone.version_id.id, ids)

    def test__get_all_user_collections_for_index_requires_auth(self):
        # test that non logged user returns 401
        test_url = furl(path="/dp/v1/user-collections/index")
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 401)

    def test__get_all_user_collections_for_index(self):
        """
        The `my-collections/index` endpoint should return all public collections and the private
        collections the user has WRITE access to. Curators can see all public collectoins and private collections where
        they are the owner. Super Curators can see all public and private collections.
        """

        public_collection = self.generate_published_collection(owner="test_user_id")
        private_collection = self.generate_unpublished_collection("test_user_id")

        public_collection_not_owned = self.generate_published_collection("test_user_id_2")
        private_collection_not_owned = self.generate_unpublished_collection("test_user_id_2")

        collection_to_tombstone = self.generate_published_collection("test_user_id")
        self.business_logic.tombstone_collection(collection_to_tombstone.collection_id)

        # test that super curators can see the all public and private collections
        test_url = furl(path="/dp/v1/user-collections/index")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token("super")}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [collection["id"] for collection in body]

        collections_by_id = {collection["id"]: collection for collection in body}

        self.assertIn(public_collection.collection_id.id, collections_by_id)
        self.assertEqual(
            collections_by_id[public_collection.collection_id.id]["curator_name"], public_collection.curator_name
        )

        self.assertIn(private_collection.collection_id.id, ids)
        self.assertEqual(
            collections_by_id[private_collection.collection_id.id]["curator_name"], private_collection.curator_name
        )

        self.assertIn(private_collection_not_owned.collection_id.id, collections_by_id)
        self.assertEqual(
            collections_by_id[private_collection_not_owned.collection_id.id]["curator_name"],
            private_collection_not_owned.curator_name,
        )

        self.assertNotIn(private_collection.version_id.id, collections_by_id)
        self.assertNotIn(collection_to_tombstone.collection_id.id, collections_by_id)
        self.assertNotIn(collection_to_tombstone.version_id.id, collections_by_id)
        self.assertNotIn(private_collection_not_owned.version_id.id, collections_by_id)

        # test that super curators can WRITE all collections
        for collection in body:
            self.assertEqual(collection["access_type"], "WRITE")

        # test that the owners can see their private collections
        # but not the private collections of other users, and that they can see public
        # collections of other users.
        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token(),
        }
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [collection["id"] for collection in body]

        self.assertIn(public_collection.collection_id.id, ids)
        self.assertIn(public_collection_not_owned.collection_id.id, ids)
        self.assertIn(private_collection.collection_id.id, ids)

        self.assertNotIn(private_collection_not_owned.collection_id.id, ids)
        self.assertNotIn(private_collection_not_owned.version_id.id, ids)
        self.assertNotIn(private_collection.version_id.id, ids)
        self.assertNotIn(collection_to_tombstone.collection_id.id, ids)
        self.assertNotIn(collection_to_tombstone.version_id.id, ids)

        # test that the owner of a collection can WRITE their collections
        # but not other users collections
        for collection in body:
            if collection["owner"] == "test_user_id":
                self.assertEqual(collection["access_type"], "WRITE")
            else:
                self.assertEqual(collection["access_type"], "READ")

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
                headers = {
                    "host": "localhost",
                    "Content-Type": "application/json",
                    "Cookie": self.get_cxguser_token(),
                }
                response = self.app.post("/dp/v1/collections", headers=headers, data=json.dumps(body))
                self.assertEqual(400, response.status_code)
                for error in expected_errors:
                    self.assertIn(error, response.json["detail"])


# 🔴 TODO: This should be reviewed. Collection deletion is a weird beast
class TestCollectionDeletion(BaseAPIPortalTest):
    def test_delete_private_collection_version__ok(self):
        # Generate test collection
        collection = self.generate_unpublished_collection()

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        test_url = furl(path=f"/dp/v1/collections/{collection.version_id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 200)

        # delete collection
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 204)

        # check collection version was deleted
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_delete_collection_revision__ok(self):
        # Generate test collection
        collection = self.generate_published_collection()
        # Generate the public collection with the same id as the private so a tombstone is created
        revision = self.generate_revision(collection.collection_id)
        dataset = self.generate_dataset(collection_version=revision)

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        test_private_url = furl(path=f"/dp/v1/collections/{revision.version_id}")
        test_public_url = furl(path=f"/dp/v1/collections/{collection.collection_id}")
        response = self.app.get(test_private_url.url, headers=headers)
        self.assertEqual(200, response.status_code)

        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset.dataset_version_id, dataset_ids)

        # delete collection revision
        response = self.app.delete(test_private_url.url, headers=headers)

        self.assertEqual(response.status_code, 204)

        # check collection revision deleted
        response = self.app.get(test_private_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

        # Public collection still exists, and does not contain revision datasets
        response = self.app.get(test_public_url.url, headers=headers)
        self.assertEqual(response.status_code, 200)

        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertNotIn(dataset.dataset_version_id, dataset_ids)

    def test_delete_collection_version__public__403(self):
        collection = self.generate_published_collection()

        test_url = furl(path=f"/dp/v1/collections/{collection.version_id}")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 405)

    def test_delete_published_collection__ok(self):
        """Published collections are tombstoned."""
        # Generate the public collection
        collection = self.generate_published_collection()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        test_public_url = furl(
            path=f"/dp/v1/collections/{collection.collection_id}", query_params=dict(visibility="PUBLIC")
        )
        # tombstone public collection
        self.business_logic.tombstone_collection(collection.collection_id)

        # check collection is gone
        response = self.app.get(test_public_url.url, headers=headers)
        self.assertEqual(response.status_code, 410)

        # check collection version also returns 'gone'
        test_version_url = furl(
            path=f"/dp/v1/collections/{collection.version_id}", query_params=dict(visibility="PUBLIC")
        )
        response = self.app.get(test_version_url.url, headers=headers)
        self.assertEqual(response.status_code, 410)
        # body = json.loads(response.data)
        # datasets_tombstoned = [dataset["tombstone"] for dataset in body["datasets"]]
        # self.assertTrue(all(datasets_tombstoned))

        # check that tombstoned collection doesn't appear in collections list endpoint
        test_list_url = furl(path="/dp/v1/collections")
        response = self.app.get(test_list_url.url, headers=headers)
        body = json.loads(response.data)
        collection_ids = [collection.id for collection in body["collections"]]
        self.assertNotIn(collection.collection_id, collection_ids)
        self.assertNotIn(collection.version_id, collection_ids)

    def test_delete_collection_version__already_deleted__403(self):
        collection = self.generate_unpublished_collection()

        # delete private collection
        test_private_url = furl(
            path=f"/dp/v1/collections/{collection.version_id}", query_params=dict(visibility="PRIVATE")
        )
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_private_url.url, headers=headers)
        self.assertEqual(response.status_code, 204)

        # check that DELETE on an already deleted collection is forbidden
        response = self.app.delete(test_private_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_delete_collection__not_owner(self):
        collection = self.generate_unpublished_collection(owner="not_test_user_id")
        test_url = furl(path=f"/dp/v1/collections/{collection.version_id}")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_delete_collection__does_not_exist(self):
        fake_id = CollectionId()
        test_url = furl(path=f"/dp/v1/collections/{fake_id}", query_params=dict(visibility="PRIVATE"))
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test_deleted_collection_does_not_appear_in_collection_lists(self):
        private_collection = self.generate_unpublished_collection()
        public_collection = self.generate_published_collection()
        collection_to_delete = self.generate_unpublished_collection()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.get("/dp/v1/collections", headers=headers)

        collection_ids = [collection["id"] for collection in json.loads(response.data)["collections"]]
        self.assertIn(private_collection.collection_id.id, collection_ids)
        self.assertIn(public_collection.collection_id.id, collection_ids)
        self.assertIn(collection_to_delete.collection_id.id, collection_ids)

        test_url = furl(
            path=f"/dp/v1/collections/{collection_to_delete.version_id.id}", query_params=dict(visibility="PRIVATE")
        )
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 204)

        # check not returned privately
        response = self.app.get("/dp/v1/collections", headers=headers)
        collection_ids = [collection["id"] for collection in json.loads(response.data)["collections"]]
        self.assertIn(private_collection.collection_id.id, collection_ids)
        self.assertIn(public_collection.collection_id.id, collection_ids)
        self.assertNotIn(collection_to_delete.version_id.id, collection_ids)
        self.assertNotIn(collection_to_delete.collection_id.id, collection_ids)

        # check not returned publicly
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get("/dp/v1/collections", headers=headers)
        collection_ids = [collection["id"] for collection in json.loads(response.data)["collections"]]
        self.assertIn(public_collection.collection_id.id, collection_ids)
        self.assertNotIn(private_collection.collection_id.id, collection_ids)
        self.assertNotIn(collection_to_delete.version_id.id, collection_ids)
        self.assertNotIn(collection_to_delete.collection_id.id, collection_ids)


class TestUpdateCollection(BaseAPIPortalTest):
    # ✅
    def test__update_collection__OK(self):
        collection = self.generate_unpublished_collection()
        test_fields = [
            "name",
            "description",
            "contact_name",
            "contact_email",
            "links",
            "consortia",
        ]
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # Update the collection
        expected_body = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
            "consortia": ["Consortia 1"],
        }
        data = json.dumps(expected_body)
        response = self.app.put(f"/dp/v1/collections/{collection.version_id.id}", data=data, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        for field in test_fields:
            self.assertEqual(expected_body[field], actual_body[field])

    def test__update_collection_strip_string_fields__OK(self):
        collection = self.generate_unpublished_collection()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # Update the collection
        new_body = {
            "name": "collection name    ",
            "description": "    This is a test collection",
            "contact_name": "   person human",
            "contact_email": "  person@human.com  ",
            "links": [{"link_name": " DOI Link ", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}],
            "consortia": ["  Consortia 1   "],
        }
        data = json.dumps(new_body)
        response = self.app.put(f"/dp/v1/collections/{collection.version_id.id}", data=data, headers=headers)

        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        self.assertEqual(new_body["name"].strip(), actual_body["name"])
        self.assertEqual(new_body["description"].strip(), actual_body["description"])
        self.assertEqual(new_body["contact_name"].strip(), actual_body["contact_name"])
        self.assertEqual(new_body["contact_email"].strip(), actual_body["contact_email"])
        self.assertEqual(["Consortia 1"], actual_body["consortia"])
        self.assertEqual(
            [{"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"}], actual_body["links"]
        )

    def test__update_collection_partial__OK(self):
        collection = self.generate_unpublished_collection(links=[Link("Link 1", "DOI", "http://doi.org/123")])
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        payload = {
            "name": "new collection name",
        }

        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id.id}", data=json.dumps(payload), headers=headers
        )
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)

        self.assertEqual(actual_body["name"], "new collection name")
        self.assertEqual(actual_body["description"], collection.metadata.description)
        self.assertEqual(actual_body["contact_name"], collection.metadata.contact_name)
        self.assertEqual(actual_body["contact_email"], collection.metadata.contact_email)
        self.assertEqual(actual_body["consortia"], collection.metadata.consortia)
        self.assertEqual(
            actual_body["links"],
            [
                {"link_name": link.name, "link_type": link.type, "link_url": link.uri}
                for link in collection.metadata.links
            ],
        )

    # ✅
    def test__update_collection__403(self):
        collection = self.generate_unpublished_collection(owner="someone else")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        data = json.dumps({"name": "new name"})
        response = self.app.put(f"/dp/v1/collections/{collection.version_id.id}", data=data, headers=headers)
        self.assertEqual(403, response.status_code)

    def test__update_collection_consortia(self):
        collection = self.generate_unpublished_collection()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # Change Consortia
        payload = {
            "consortia": ["Consortia 1", "Consortia 3"],
        }

        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id.id}", data=json.dumps(payload), headers=headers
        )
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)

        self.assertEqual(actual_body["name"], collection.metadata.name)
        self.assertEqual(actual_body["description"], collection.metadata.description)
        self.assertEqual(actual_body["contact_name"], collection.metadata.contact_name)
        self.assertEqual(actual_body["contact_email"], collection.metadata.contact_email)
        self.assertEqual(actual_body["consortia"], ["Consortia 1", "Consortia 3"])

    def test__remove_collection_consortia(self):
        collection = self.generate_unpublished_collection()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # Remove Consortia
        payload = {
            "consortia": ["Consortia 1"],
        }

        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id.id}", data=json.dumps(payload), headers=headers
        )
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)

        self.assertEqual(actual_body["name"], collection.metadata.name)
        self.assertEqual(actual_body["description"], collection.metadata.description)
        self.assertEqual(actual_body["contact_name"], collection.metadata.contact_name)
        self.assertEqual(actual_body["contact_email"], collection.metadata.contact_email)
        self.assertEqual(actual_body["consortia"], ["Consortia 1"])

    def test__remove_all_collection_consortia(self):
        collection = self.generate_unpublished_collection()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # Remove All Consortia
        payload = {
            "consortia": [],
        }

        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id.id}", data=json.dumps(payload), headers=headers
        )
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)

        self.assertEqual(actual_body["name"], collection.metadata.name)
        self.assertEqual(actual_body["description"], collection.metadata.description)
        self.assertEqual(actual_body["contact_name"], collection.metadata.contact_name)
        self.assertEqual(actual_body["contact_email"], collection.metadata.contact_email)
        self.assertEqual(actual_body["consortia"], [])

    def test__add_new_collection_consortia(self):
        collection = self.generate_unpublished_collection()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # Add Consortia
        payload = {
            "consortia": ["Consortia 4"],
        }

        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id.id}", data=json.dumps(payload), headers=headers
        )
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)

        self.assertEqual(actual_body["name"], collection.metadata.name)
        self.assertEqual(actual_body["description"], collection.metadata.description)
        self.assertEqual(actual_body["contact_name"], collection.metadata.contact_name)
        self.assertEqual(actual_body["contact_email"], collection.metadata.contact_email)
        self.assertEqual(actual_body["consortia"], ["Consortia 4"])

    def test__update_with_invalid_collection_consortia(self):
        collection = self.generate_unpublished_collection()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # Invalid Consortia
        payload = {
            "consortia": ["Invalid Consortia"],
        }

        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id.id}", data=json.dumps(payload), headers=headers
        )
        self.assertEqual(400, response.status_code)
        error_payload = json.loads(response.data)
        self.assertEqual([{"name": "consortia", "reason": "Invalid consortia."}], error_payload["detail"])

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
        collection = self.generate_unpublished_collection(links=[Link("Link 1", "DOI", "http://doi.org/123")])

        self.assertIsNotNone(collection.publisher_metadata)
        if collection.publisher_metadata:  # pylance
            self.assertEqual("Old Journal", collection.publisher_metadata["journal"])

        # From now on, Crossref will return `New Journal`
        self.crossref_provider.fetch_metadata = Mock(return_value=generate_mock_publisher_metadata("New Journal"))

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.put(
            f"/dp/v1/collections/{collection.version_id}",
            data=json.dumps({"links": [{"link_name": "Link 1", "link_url": "10.1234/5678", "link_type": "DOI"}]}),
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
        collection = self.generate_unpublished_collection(links=[Link("Link 1", "DOI", "http://doi.org/123")])

        self.assertIsNotNone(collection.publisher_metadata)
        if collection.publisher_metadata:  # pylance
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

        collection = self.generate_unpublished_collection(links=[Link("Link 1", "DOI", "http://doi.org/123")])

        self.assertIsNotNone(collection.publisher_metadata)
        if collection.publisher_metadata:  # pylance
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


class TestCollectionsCurators(BaseAPIPortalTest):
    def test_view_non_owned_private_collection__ok(self):
        # Generate test collection
        collection = self.generate_unpublished_collection(owner="another_test_user_id")

        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token(user="not_owner"),
        }
        test_url = furl(path=f"/dp/v1/collections/{collection.version_id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers=headers)

        # This will pass even for non curators.
        # Why are users allowed to view private collections that they don't own?
        self.assertEqual(response.status_code, 200)

        body = json.loads(response.data)
        self.assertEqual(body["access_type"], "READ")

    # ✅
    def test_update_non_owned_private_collection_as_super_curator__ok(self):
        # Generate test collection
        collection = self.generate_unpublished_collection(owner="another_test_user_id")

        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token("super"),
        }

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

    def test_delete_non_owned_private_collection_as_super_curator__ok(self):
        # Generate test collection
        collection = self.generate_unpublished_collection(owner="another_test_user_id")

        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token("super"),
        }

        test_url = furl(path=f"/dp/v1/collections/{collection.version_id}", query_params=dict(visibility="PRIVATE"))
        response = self.app.delete(test_url.url, headers=headers)
        self.assertEqual(204, response.status_code)


# TODO: these tests all require the generation of a dataset
class TestDataset(BaseAPIPortalTest):
    # ✅
    def test__post_dataset_asset__OK(self):
        self.business_logic.get_dataset_artifact_download_data = Mock(
            return_value=DatasetArtifactDownloadData(
                "asset.h5ad", DatasetArtifactType.H5AD, 1000, "http://presigned.url"
            )
        )
        version = self.generate_dataset(
            artifacts=[DatasetArtifactUpdate(DatasetArtifactType.H5AD, "http://mock.uri/asset.h5ad")]
        )
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
        version = self.generate_dataset(
            artifacts=[DatasetArtifactUpdate(DatasetArtifactType.H5AD, "http://mock.uri/asset.h5ad")]
        )
        dataset_version_id = version.dataset_version_id
        artifact_id = version.artifact_ids[0]

        test_url = furl(path=f"/dp/v1/datasets/{dataset_version_id}/asset/{artifact_id}")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(500, response.status_code)

    # ✅
    def test__post_dataset_asset__dataset_NOT_FOUND(self):
        fake_id = DatasetVersionId()
        test_url = furl(path=f"/dp/v1/datasets/{fake_id}/asset/{fake_id}")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.data)
        self.assertEqual(f"'dataset/{fake_id}' not found.", body["detail"])

    # ✅
    def test__post_dataset_asset__asset_NOT_FOUND(self):
        dataset = self.generate_dataset()
        test_url = furl(path=f"/dp/v1/datasets/{dataset.dataset_version_id}/asset/fake_asset")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.data)
        self.assertEqual(f"'dataset/{dataset.dataset_version_id}/asset/fake_asset' not found.", body["detail"])

    # ✅
    def test__get_status__ok(self):
        dataset = self.generate_dataset(
            statuses=[
                DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.PENDING),
                DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING),
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
            "id": "NA",  # TODO: I am deprecating this, I don't think it has any use.
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
        # TODO: decide whether to use dataset_id or version_id in this endpoint

        private_dataset = self.generate_dataset()
        public_dataset = self.generate_dataset(publish=True)

        test_url = furl(path="/dp/v1/datasets/index")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [d["id"] for d in body]
        self.assertIn(public_dataset.dataset_version_id, ids)
        self.assertNotIn(private_dataset.dataset_version_id, ids)

        actual_dataset = None
        for d in body:
            if d["id"] == public_dataset.dataset_version_id:
                actual_dataset = d

        persisted_dataset = self.business_logic.get_dataset_version(DatasetVersionId(public_dataset.dataset_version_id))
        self.assertIsNotNone(actual_dataset)
        self.assertIsNotNone(persisted_dataset)

        def convert_ontology(ontologies):
            return [dataclasses.asdict(o) for o in ontologies]

        if actual_dataset is not None and persisted_dataset is not None:  # pylance
            self.assertNotIn("description", actual_dataset)
            self.assertEqual(actual_dataset["id"], persisted_dataset.version_id.id)
            self.assertEqual(actual_dataset["name"], persisted_dataset.metadata.name)
            # self.assertEqual(actual_dataset["collection_id"], persisted_dataset.collection_id)
            self.assertEqual(actual_dataset["assay"], convert_ontology(persisted_dataset.metadata.assay))
            self.assertEqual(actual_dataset["tissue"], convert_ontology(persisted_dataset.metadata.tissue))
            self.assertEqual(actual_dataset["disease"], convert_ontology(persisted_dataset.metadata.disease))
            self.assertEqual(actual_dataset["sex"], convert_ontology(persisted_dataset.metadata.sex))
            self.assertEqual(
                actual_dataset["self_reported_ethnicity"],
                convert_ontology(persisted_dataset.metadata.self_reported_ethnicity),
            )
            self.assertEqual(actual_dataset["organism"], convert_ontology(persisted_dataset.metadata.organism))
            self.assertEqual(
                actual_dataset["development_stage"], convert_ontology(persisted_dataset.metadata.development_stage)
            )
            self.assertEqual(actual_dataset["cell_count"], persisted_dataset.metadata.cell_count)
            self.assertEqual(actual_dataset["primary_cell_count"], persisted_dataset.metadata.primary_cell_count)
            self.assertEqual(actual_dataset["cell_type"], convert_ontology(persisted_dataset.metadata.cell_type))
            self.assertEqual(actual_dataset["is_primary_data"], persisted_dataset.metadata.is_primary_data)
            self.assertEqual(actual_dataset["mean_genes_per_cell"], persisted_dataset.metadata.mean_genes_per_cell)
            # self.assertEqual(actual_dataset["explorer_url"], persisted_dataset.explorer_url)
            self.assertEqual(
                actual_dataset["published_at"], persisted_dataset.canonical_dataset.published_at.timestamp()
            )
            # self.assertEqual(actual_dataset["revised_at"], persisted_dataset.revised_at.timestamp())

    def test__get_all_datasets_for_user_index(self):
        # return the datasets related to the collections the user has WRITE access to.

        private_dataset = self.generate_dataset()
        public_dataset = self.generate_dataset(publish=True)
        self.business_logic.create_collection_version(CollectionId(public_dataset.collection_id))

        test_url = furl(path="/dp/v1/user-datasets/index")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        print(body)

        # initialize an empty list to store the dataset ids
        dataset_ids = []

        # iterate over the collections in the body list
        for collection in body:
            # iterate over the datasets in the "dataset_assets" list for each collection
            for dataset in collection["dataset_assets"]:
                # append the dataset_id to the list of dataset ids
                dataset_ids.append(dataset["dataset_id"])

        self.assertEqual(len(dataset_ids), 6)

        dataset_ids = list(dataset_ids)
        print("dataset_ids: ", dataset_ids)
        print("public_dataset.dataset_version_id: ", public_dataset.dataset_version_id)
        print("private_dataset.dataset_version_id: ", private_dataset.dataset_version_id)

        self.assertIn(public_dataset.dataset_version_id, dataset_ids)
        self.assertIn(private_dataset.dataset_version_id, dataset_ids)

    def test__get_all_user_datasets_for_index_requires_auth(self):
        # test that non logged user returns 401
        test_url = furl(path="/dp/v1/user-datasets/index")
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(response.status_code, 401)

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
            if d["id"] == dataset.dataset_version_id:
                actual_dataset = d
        self.assertIsNotNone(actual_dataset)

        def convert_ontology(ontologies):
            return [dataclasses.asdict(o) for o in ontologies]

        if actual_dataset is not None:  # pylance
            self.assertEqual(actual_dataset["development_stage"], convert_ontology(modified_metadata.development_stage))
            self.assertEqual(
                actual_dataset["development_stage_ancestors"],
                ["HsapDv:0000008", "HsapDv:0000006", "HsapDv:0000002", "HsapDv:0000045", "HsapDv:0000001"],
            )

            self.assertEqual(actual_dataset["tissue"], convert_ontology(modified_metadata.tissue))
            self.assertCountEqual(
                actual_dataset["tissue_ancestors"],
                [
                    "UBERON:0005178",
                    "UBERON:0000072",
                    "UBERON:0001558",
                    "UBERON:0000915",
                    "UBERON:0001005",
                    "UBERON:0005181",
                    "UBERON:0002075",
                    "UBERON:0000170",
                    "UBERON:0002048",
                    "UBERON:0002100",
                    "UBERON:0000171",
                    "UBERON:0009569",
                    "UBERON:0000065",
                    "UBERON:0005177",
                    "UBERON:0001004",
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
        dataset = self.generate_dataset(
            artifacts=[
                DatasetArtifactUpdate(DatasetArtifactType.CXG, "s3://mock-bucket/mock-key.cxg"),
                DatasetArtifactUpdate(DatasetArtifactType.H5AD, "s3://mock-bucket/mock-key.h5ad"),
            ]
        )

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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.WAITING)]
        )
        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # Test while uploading
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)]
        )
        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

    # ✅
    def test__cancel_dataset_download__dataset_does_not_exist(self):
        fake_id = DatasetVersionId()
        test_url = f"/dp/v1/datasets/{fake_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    # ✅
    def test__delete_uploaded_dataset__ok(self):
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)]
        )
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}

        # check dataset in collection
        collection_url = furl(path=f"/dp/v1/collections/{dataset.collection_version_id}")
        # TODO: do we still have a query for private/public?
        response = self.app.get(collection_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset.dataset_version_id, dataset_ids)

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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)]
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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADED)],
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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.WAITING)],
        )

        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    # ✅
    def test__cancel_dataset_download__user_not_logged_in(self):
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.WAITING)],
        )

        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 401)

    # ✅
    def test__dataset_meta__ok(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}

        # with self.subTest("dataset is public"):
        #     test_uri_0 = "s3://some_bucket/some_key_0.cxg"

        #     public_dataset = self.generate_dataset(
        #         artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri_0)],
        #         publish=True,
        #     )

        #     test_url_public = f"/dp/v1/datasets/meta?url={public_dataset.explorer_url}"
        #     response = self.app.get(test_url_public, headers)
        #     self.assertEqual(response.status_code, 200)

        #     expected_identifiers = {
        #         "s3_uri": test_uri_0,
        #         "dataset_id": public_dataset.dataset_version_id,
        #         "collection_id": public_dataset.collection_version_id,
        #         "collection_visibility": "PUBLIC",  # this is a published collection
        #         "tombstoned": False,
        #     }

        #     self.assertEqual(json.loads(response.data), expected_identifiers)

        # with self.subTest("dataset is private"):
        #     test_uri_1 = "s3://some_bucket/some_key_1.cxg"

        #     public_dataset = self.generate_dataset(
        #         artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri_1)],
        #         publish=False,
        #     )

        #     test_url_private = f"/dp/v1/datasets/meta?url={public_dataset.explorer_url}"
        #     expected_identifiers = {
        #         "s3_uri": test_uri_1,
        #         "dataset_id": public_dataset.dataset_version_id,
        #         "collection_id": public_dataset.collection_version_id,
        #         "collection_visibility": "PRIVATE",
        #         "tombstoned": False,
        #     }

        #     response = self.app.get(test_url_private, headers)
        #     self.assertEqual(response.status_code, 200)

        #     self.assertEqual(json.loads(response.data), expected_identifiers)

        with self.subTest("dataset is withdrawn"):
            test_uri_1 = "s3://some_bucket/some_key_1.cxg"

            public_dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri_1)],
                publish=True,
            )

            self.business_logic.tombstone_collection(CollectionId(public_dataset.collection_id))

            test_url_private = f"/dp/v1/datasets/meta?url={public_dataset.explorer_url}"
            expected_identifiers = {
                "s3_uri": test_uri_1,
                "dataset_id": public_dataset.dataset_version_id,
                "collection_id": public_dataset.collection_version_id,
                "collection_visibility": "PUBLIC",
                "tombstoned": True,
            }

            response = self.app.get(test_url_private, headers)
            self.assertEqual(response.status_code, 200)

            self.assertEqual(json.loads(response.data), expected_identifiers)

    # ✅
    def test__dataset_meta__404(self):
        fake_id = DatasetVersionId()
        headers = {"host": "localhost", "Content-Type": "application/json"}
        fake_url = f"http://base.url/{fake_id}.cxg/"
        test_url_404 = f"/dp/v1/datasets/meta?url={fake_url}"

        response = self.app.get(test_url_404, headers)
        self.assertEqual(response.status_code, 404)

    def test__explorer_portal_integration(self):
        """
        Tests the explorer <-> portal integration.
        The steps carried out by this test are:
        1. Generate the explorer_url
        2. Call the `get_dataset_identifiers` endpoint, retrieve `collection_id` and `dataset_id` from there
        3. Call the GET /collections/:collection_id endpoint, locate the dataset
        """
        headers = {"host": "localhost", "Content-Type": "application/json"}

        def _call_meta_endpoint(explorer_url):
            test_url = f"/dp/v1/datasets/meta?url={explorer_url}"
            response = self.app.get(test_url, headers)
            self.assertEqual(response.status_code, 200)
            return json.loads(response.data)

        def _call_collections_endpoint(collection_id):
            test_url = f"/dp/v1/collections/{collection_id}"
            response = self.app.get(test_url, headers)
            self.assertEqual(response.status_code, 200)
            return json.loads(response.data)

        with self.subTest("Dataset belonging to an unpublished collection"):
            test_uri = "s3://some_bucket/some_key_0.cxg"

            dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri)],
                publish=False,
            )
            # In this case, explorer_url points to the canonical link
            explorer_url = f"http://base.url/{dataset.dataset_id}.cxg/"
            meta_response = _call_meta_endpoint(explorer_url)

            returned_collection_id = meta_response["collection_id"]
            returned_dataset_id = meta_response["dataset_id"]

            collections_response = _call_collections_endpoint(returned_collection_id)
            datasets = collections_response["datasets"]
            self.assertIn(returned_dataset_id, [dataset["id"] for dataset in datasets])

        with self.subTest("Dataset belonging to an unpublished collection, replaced"):
            test_uri = "s3://some_bucket/some_key_0.cxg"

            dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri)],
                publish=False,
            )

            collection_version = self.business_logic.get_collection_version(
                CollectionVersionId(dataset.collection_version_id)
            )

            replaced_dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri)],
                collection_version=collection_version,
                replace_dataset_version_id=DatasetVersionId(dataset.dataset_version_id),
            )
            # In this case, explorer_url points to the canonical link
            explorer_url = f"http://base.url/{dataset.dataset_id}.cxg/"
            meta_response = _call_meta_endpoint(explorer_url)

            returned_collection_id = meta_response["collection_id"]
            returned_dataset_id = meta_response["dataset_id"]

            self.assertEqual(returned_dataset_id, replaced_dataset.dataset_version_id)

            collections_response = _call_collections_endpoint(returned_collection_id)
            datasets = collections_response["datasets"]
            self.assertIn(returned_dataset_id, [dataset["id"] for dataset in datasets])

        with self.subTest("Dataset belonging to a published collection"):
            test_uri = "s3://some_bucket/some_key_1.cxg"

            dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri)], publish=True
            )
            # In this case, explorer_url points to the canonical link
            explorer_url = f"http://base.url/{dataset.dataset_id}.cxg/"
            meta_response = _call_meta_endpoint(explorer_url)

            returned_collection_id = meta_response["collection_id"]
            returned_dataset_id = meta_response["dataset_id"]

            collections_response = _call_collections_endpoint(returned_collection_id)
            datasets = collections_response["datasets"]
            self.assertIn(returned_dataset_id, [dataset["id"] for dataset in datasets])

        with self.subTest("Dataset belonging to a revision of a published collection, not replaced"):
            test_uri = "s3://some_bucket/some_key_2.cxg"

            dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri)], publish=True
            )
            self.business_logic.create_collection_version(CollectionId(dataset.collection_id))

            # In this case, explorer_url points to the versioned link
            explorer_url = f"http://base.url/{dataset.dataset_version_id}.cxg/"
            meta_response = _call_meta_endpoint(explorer_url)

            returned_collection_id = meta_response["collection_id"]
            returned_dataset_id = meta_response["dataset_id"]

            collections_response = _call_collections_endpoint(returned_collection_id)
            datasets = collections_response["datasets"]
            self.assertIn(returned_dataset_id, [dataset["id"] for dataset in datasets])

        with self.subTest("Dataset belonging to a revision of a published collection, replaced"):
            test_uri = "s3://some_bucket/some_key_1.cxg"

            dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri)], publish=True
            )
            revision = self.business_logic.create_collection_version(CollectionId(dataset.collection_id))
            revised_dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri)],
                collection_version=revision,
                replace_dataset_version_id=DatasetVersionId(dataset.dataset_version_id),
            )
            self.assertEqual(revised_dataset.dataset_id, dataset.dataset_id)
            self.assertNotEqual(revised_dataset.dataset_version_id, dataset.dataset_version_id)

            # Retrieve the explorer url from the GET collections/:collection_id endpoint. This is the only way to force
            # explorer_url to be exactly the same used by the portal to open the explorer url
            test_url = f"/dp/v1/collections/{revision.version_id}"
            response = self.app.get(test_url, headers)
            self.assertEqual(response.status_code, 200)
            response_data = json.loads(response.data)
            datasets = response_data["datasets"]
            self.assertIn(revised_dataset.dataset_version_id, [dataset["id"] for dataset in datasets])
            replaced_dataset = next(
                dataset for dataset in datasets if dataset["id"] == revised_dataset.dataset_version_id
            )

            explorer_url = replaced_dataset["dataset_deployments"][0]["url"]
            meta_response = _call_meta_endpoint(explorer_url)

            returned_collection_id = meta_response["collection_id"]
            returned_dataset_id = meta_response["dataset_id"]

            collections_response = _call_collections_endpoint(returned_collection_id)
            datasets = collections_response["datasets"]
            self.assertIn(returned_dataset_id, [dataset["id"] for dataset in datasets])

        with self.subTest("Dataset that appears in multiple published versions"):
            """
            If a dataset appears in multiple collection versions, the most recent one will be returned.
            """
            test_uri = "s3://some_bucket/some_key_1.cxg"

            dataset = self.generate_dataset(
                artifacts=[DatasetArtifactUpdate(DatasetArtifactType.CXG, test_uri)], publish=True
            )
            revision = self.business_logic.create_collection_version(CollectionId(dataset.collection_id))

            self.business_logic.publish_collection_version(revision.version_id)

            # Both versions are now published
            original_version = self.business_logic.get_collection_version(
                CollectionVersionId(dataset.collection_version_id)
            )
            revision_version = self.business_logic.get_collection_version(revision.version_id)

            self.assertIsNotNone(original_version.published_at)
            self.assertIsNotNone(revision_version.published_at)

            explorer_url = f"http://base.url/{dataset.dataset_version_id}.cxg/"
            meta_response = _call_meta_endpoint(explorer_url)

            returned_collection_id = meta_response["collection_id"]
            returned_dataset_id = meta_response["dataset_id"]

            self.assertEqual(returned_collection_id, revision_version.version_id.id)

            collections_response = _call_collections_endpoint(returned_collection_id)
            datasets = collections_response["datasets"]
            self.assertIn(returned_dataset_id, [dataset["id"] for dataset in datasets])

        with self.subTest("Dataset that is not found"):
            explorer_url = "http://base.url/01234567-89ab-cdef-0123-456789abcdef.cxg/"
            test_url = f"/dp/v1/datasets/meta?url={explorer_url}"
            res = self.app.get(test_url, headers)
            self.assertEqual(404, res.status_code)


class TestDatasetCurators(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()

    def tearDown(self):
        super().tearDown()

    def test__get_status__200_for_non_owned_dataset_as_super_curator(self):
        processing_status = DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.WAITING)
        dataset = self.generate_dataset(owner="someone_else", statuses=[processing_status])
        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}/status"
        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token("super"),
        }
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(200, response.status_code)

    # 💛 review later
    def test__cancel_dataset_download__202_user_not_collection_owner_as_super_curator(self):
        processing_status = DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.WAITING)
        dataset = self.generate_dataset(owner="someone_else", statuses=[processing_status])
        test_url = f"/dp/v1/datasets/{dataset.dataset_version_id}"
        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token("super"),
        }
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)


# #### REVISIONS START HERE #####


class TestRevision(BaseAPIPortalTest):
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
        self.assertCountEqual([link["link_name"] for link in response_post_json["links"]], ["Link 1", "DOI Link"])

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
        fake_collection_id = CollectionId()
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        response = self.app.post(f"/dp/v1/collections/{fake_collection_id}", headers=headers)
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


class TestDeleteRevision(BaseAPIPortalTest):
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
        published_collection = self.generate_published_collection(owner="someone_else", add_datasets=2)
        revision_id = self._create_revision(published_collection.collection_id.id, user="super")

        # Delete the revision
        path = f"/dp/v1/collections/{revision_id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        resp = self.app.delete(path, headers=headers)
        self.assertEqual(403, resp.status_code)


# Those were the previous revision tests - mostly covered by the business logic layer now
# A few cases are good, but we don't need to test every single case here
class TestPublishRevision(BaseAPIPortalTest):
    """Test case for publishing a revision."""

    def setUp(self):
        super().setUp()
        self.base_path = "/dp/v1/collections"
        self.mock_timestamp = datetime(2000, 12, 25, 0, 0)
        # TODO: header pattern
        self.headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token(),
        }

    # ✅
    def test__publish_revision_with_new_dataset__OK(self):
        """
        Publish a revision with new datasets.
        """

        unpublished_collection = self.generate_unpublished_collection(add_datasets=1)

        path = f"{self.base_path}/{unpublished_collection.version_id.id}/publish"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
        body = {"data_submission_policy_version": "1.0"}  # TODO: still in use?
        response = self.app.post(path, headers=headers, data=json.dumps(body))

        self.assertEqual(202, response.status_code)
        self.assertDictEqual(
            {"collection_id": unpublished_collection.collection_id.id, "visibility": "PUBLIC"},
            json.loads(response.data),
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

    # ✅
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
        body = {"data_submission_policy_version": "1.0"}  # TODO: still in use?

        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)

        response_json = json.loads(response.data)
        collection_id = response_json["collection_id"]

        # Checks that the published revision exists and only has one dataset
        published_collection = self.business_logic.get_published_collection_version(CollectionId(collection_id))
        self.assertIsNotNone(published_collection)
        self.assertIsNotNone(published_collection.published_at)
        self.assertEqual(published_collection.collection_id, unpublished_collection.collection_id)
        self.assertEqual(1, len(published_collection.datasets))

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
        body = {"data_submission_policy_version": "1.0"}  # TODO: still in use?
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(409, response.status_code)

    def test__publish_revision_as_non_autorized__401(self):
        """
        Publish a collection as a non authenticated user returns 403
        """
        collection = self.generate_unpublished_collection(add_datasets=1)
        path = f"{self.base_path}/{collection.version_id.id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(401, response.status_code)

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
        collection_id = CollectionId()
        body = {"data_submission_policy_version": "1.0"}
        path = f"{self.base_path}/{collection_id}/publish"
        response = self.app.post(path, headers=self.headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    # ✅
    def test__can_publish_owned_collection(self):
        collection = self.generate_unpublished_collection(add_datasets=1)
        path = f"{self.base_path}/{collection.version_id.id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token("owner"),
        }
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)

    # ✅
    def test__can_publish_non_owned_collection_as_super_curator(self):
        collection = self.generate_unpublished_collection(add_datasets=1, owner="someone else")
        path = f"{self.base_path}/{collection.version_id.id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token("super"),
        }
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)

    def test__publish_collection_triggers_cloudfront_invalidation(self):
        self.cloudfront_provider.create_invalidation_for_index_paths = Mock()
        collection = self.generate_unpublished_collection(add_datasets=1)
        path = f"{self.base_path}/{collection.version_id.id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token(),
        }
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)
        self.cloudfront_provider.create_invalidation_for_index_paths.assert_called_once()

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
    #     new_dataset_id = self.generate_dataset_with_s3_resources(self.session,
    #     collection_id=self.rev_collection.id).id
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
    #     new_dataset_id = self.generate_dataset_with_s3_resources(self.session,
    #     collection_id=self.rev_collection.id).id
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


# ##### UPLOAD TESTS START HERE ######


class TestCollectionPostUploadLink(BaseAPIPortalTest):
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

        with self.subTest("Bad Dropbox link"):
            self.uri_provider.validate = Mock(return_value=True)
            self.uri_provider.get_file_info.side_effect = FileInfoException(
                "The URL provided causes an error with Dropbox."
            )
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": self.get_cxguser_token()}
            body = {"url": self.dummy_link}
            test_url = furl(path=path)
            response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
            self.assertEqual(400, response.status_code)
            print(json.loads(response.data)["detail"])
            self.assertTrue(json.loads(response.data)["detail"] == "The URL provided causes an error with Dropbox.")

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
        fake_id = CollectionId()
        path = f"/dp/v1/collections/{fake_id}/upload-links"
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


class TestCollectionPutUploadLink(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        # TODO: headers do not follow the same pattern as the rest of the test. Maybe change?
        self.headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token(),
        }
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "file.h5ad"))

    def assert_datasets_are_updated(self, dataset1, dataset2):
        """
        Asserts that dataset2 is an updated revision of dataset1
        """
        self.assertIsNotNone(dataset2)
        if dataset2 is not None:  # pylance
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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)], publish=True
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
    @patch("backend.layers.thirdparty.step_function_provider.StepFunctionProvider")
    def test__reupload_unpublished_dataset__202(self, mock_upload_sfn):
        """
        Reuploads an unpublished dataset
        """
        dataset = self.generate_dataset(
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.PENDING)],
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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
        )
        collection = self.generate_unpublished_collection()

        path = f"/dp/v1/collections/{collection.version_id.id}/upload-links"
        body = {"url": self.good_link, "id": dataset.dataset_version_id}

        response = self.app.put(path, headers=self.headers, data=json.dumps(body))
        self.assertEqual(404, response.status_code)


class TestCollectionUploadLinkCurators(BaseAPIPortalTest):
    def setUp(self):
        super().setUp()
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        self.uri_provider.get_file_info = Mock(return_value=FileInfo(1, "file.h5ad"))

    # ✅
    def test__can_upload_dataset_to_non_owned_collection_as_super_curator(self):
        """
        A super curator can upload a dataset to a non-owned collection
        """
        collection = self.generate_unpublished_collection(owner="someone else")
        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token("super"),
        }
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
            statuses=[DatasetStatusUpdate(DatasetStatusKey.PROCESSING, DatasetProcessingStatus.SUCCESS)],
        )

        path = f"/dp/v1/collections/{dataset.collection_version_id}/upload-links"
        body = {"url": self.good_link, "id": dataset.dataset_version_id}

        headers = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": self.get_cxguser_token("super"),
        }
        response = self.app.put(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)
