import itertools
import json
import sys
import os
import unittest
from datetime import datetime

from furl import furl
from multiprocessing import Process

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Collection
from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.utils import BogusCollectionParams
from tests.unit.backend.chalice.api_server.mock_auth import PORT, launch_mock_oauth, get_auth_token


class TestCollection(BaseAPITest, unittest.TestCase):
    def setUp(self):
        self.mock_oauth_process = Process(target=launch_mock_oauth)
        self.mock_oauth_process.start()

        old_path = sys.path.copy()

        def restore_path(p):
            sys.path = p

        sys.path.insert(0, os.path.join(self.corpora_api_dir, "chalicelib"))  # noqa
        self.addCleanup(restore_path, old_path)
        from corpora.common.corpora_config import CorporaAuthConfig

        # Use the CorporaAuthConfig used by the chalice app
        self.auth_config = CorporaAuthConfig()
        self.auth_config._config["api_base_url"] = f"http://localhost:{PORT}"
        self.auth_config._config["callback_base_url"] = "http://localhost:5000"
        self.auth_config.update_defaults()

    def tearDown(self):
        self.mock_oauth_process.terminate()

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

    def generate_collection(self, **params):
        _collection = Collection.create(**BogusCollectionParams.get(**params))
        # Cleanup collection after test
        self.addCleanup(_collection.delete)
        return _collection.id

    def test__list_collection__ok(self):
        path = "/dp/v1/collections"
        headers = dict(host="localhost")

        from_date = int(datetime.fromtimestamp(60).timestamp())
        creation_time = 70
        to_date = int(datetime.fromtimestamp(80).timestamp())

        expected_id = self.generate_collection(
            visibility=CollectionVisibility.PUBLIC.name, created_at=datetime.fromtimestamp(creation_time)
        )

        with self.subTest("No Parameters"):
            test_url = furl(path=path)
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.validate_collections_response_structure(actual_body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["collections"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(None, actual_body.get("from_date"))

        with self.subTest("from_date"):
            test_url = furl(path=path, query_params={"from_date": from_date})
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.validate_collections_response_structure(actual_body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["collections"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(actual_body["from_date"], from_date)

        with self.subTest("to_date"):
            test_url = furl(path=path, query_params={"to_date": to_date})
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.validate_collections_response_structure(actual_body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["collections"]])
            self.assertEqual(to_date, actual_body["to_date"])
            self.assertEqual(None, actual_body.get("from_date"))

        with self.subTest("from_date->to_date"):
            test_url = furl(path=path, query_params={"from_date": from_date, "to_date": to_date})
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
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
                    "dataset_deployments": [
                        {
                            "dataset_id": "test_dataset_id",
                            "id": "test_deployment_directory_id",
                            "url": "test_url",
                        }
                    ],
                    "development_stage": [{"label": "test_develeopment_stage", "ontology_term_id": "test_obo"}],
                    "disease": [
                        {"label": "test_disease", "ontology_term_id": "test_obo"},
                        {"label": "test_disease2", "ontology_term_id": "test_obp"},
                        {"label": "test_disease3", "ontology_term_id": "test_obq"},
                    ],
                    "ethnicity": [{"label": "test_ethnicity", "ontology_term_id": "test_obo"}],
                    "id": "test_dataset_id",
                    "name": "test_dataset_name",
                    "organism": {"label": "test_organism", "ontology_term_id": "test_obo"},
                    "collection_id": "test_collection_id",
                    "collection_visibility": "PUBLIC",
                    "cell_count": None,
                    "is_valid": False,
                    "revision": 0,
                    "sex": ["test_sex", "test_sex2"],
                    "tissue": [{"label": "test_tissue", "ontology_term_id": "test_obo"}],
                    "processing_status": {
                        "id": "test_dataset_processing_status_id",
                        "dataset_id": "test_dataset_id",
                        "upload_status": "UPLOADING",
                        "upload_message": None,
                        "upload_progress": 4 / 9,
                        "validation_status": "NA",
                        "validation_message": None,
                        "conversion_loom_status": "NA",
                        "conversion_anndata_status": "NA",
                        "conversion_rds_status": "NA",
                        "conversion_cxg_status": "NA",
                    },
                }
            ],
            "description": "test_description",
            "id": "test_collection_id",
            "links": [
                {"link_type": "RAW_DATA", "link_name": "test_link_name", "link_url": "test_url"},
                {"link_type": "OTHER", "link_name": "test_summary_link_name", "link_url": "test_summary_url"},
            ],
            "name": "test_collection",
            "visibility": "PUBLIC",
            "obfuscated_uuid": "",
            "contact_email": "",
            "contact_name": "",
            "data_submission_policy_version": "0",
        }

        with self.subTest("auth cookie"):
            expected_body["access_type"] = "WRITE"
            test_url = furl(path="/dp/v1/collections/test_collection_id", query_params=dict(visibility="PUBLIC"))
            cxguser_cookie = get_auth_token(self.app)
            response = self.app.get(test_url.url, headers=dict(host="localhost", Cookie=cxguser_cookie))
            response.raise_for_status()
            self.validate_collection_uuid_response_structure(json.loads(response.body))
            actual_body = self.remove_timestamps(json.loads(response.body))
            self.assertDictEqual(actual_body, expected_body)

        with self.subTest("no auth cookie"):
            expected_body["access_type"] = "READ"
            test_url = furl(path="/dp/v1/collections/test_collection_id", query_params=dict(visibility="PUBLIC"))
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            self.validate_collection_uuid_response_structure(json.loads(response.body))
            actual_body = self.remove_timestamps(json.loads(response.body))
            self.assertDictEqual(actual_body, expected_body)

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
                visibility=CollectionVisibility.PUBLIC.name, owner="someone else"
            ),
            public_owned=self.generate_collection(visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"),
            private_not_owner=self.generate_collection(
                visibility=CollectionVisibility.PRIVATE.name, owner="someone else"
            ),
            private_owned=self.generate_collection(visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"),
        )

        # run
        for auth, owns, visi, obfu in test_cases:
            if obfu or (visi == "private" and owns and auth) or (visi == "public"):
                expected_response_code = 200
            else:
                expected_response_code = 403

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
                    actual_body = json.loads(response.body)
                    self.assertEqual(expected_access_type, actual_body["access_type"])

    def test__get_collection_uuid__403_not_found(self):
        """Verify the test collection exists and the expected fields exist."""
        test_url = furl(path="/dp/v1/collections/AAAA-BBBB-CCCC-DDDD", query_params=dict(visibility="PUBLIC"))
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, response.status_code)
        self.assertIn("X-AWS-REQUEST-ID", response.headers.keys())

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
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
            ],
        }
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(data))
        self.assertEqual(201, response.status_code)
        collection_uuid = json.loads(response.body)["collection_uuid"]

        test_url = furl(path=f"/dp/v1/collections/{collection_uuid}")
        test_url.add(query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.body)

        self.assertEqual(body["description"], data["description"])
        self.assertEqual(body["name"], data["name"])
        self.assertEqual(body["contact_name"], body["contact_name"])
        self.assertEqual(body["contact_email"], body["contact_email"])

        # test that non owners cant access
        no_cookie_headers = {"host": "localhost", "Content-Type": "application/json"}
        test_url = furl(path=f"/dp/v1/collections/{collection_uuid}")
        test_url.add(query_params=dict(visibility="PRIVATE"))
        response = self.app.get(test_url.url, no_cookie_headers)
        self.assertEqual(403, response.status_code)

    def test__list_collection__check_owner(self):

        # Generate test collection
        public_owned = self.generate_collection(visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id")
        private_owned = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id")
        public_not_owned = self.generate_collection(visibility=CollectionVisibility.PUBLIC.name, owner="someone else")
        private_not_owned = self.generate_collection(visibility=CollectionVisibility.PRIVATE.name, owner="someone else")

        path = "/dp/v1/collections"
        with self.subTest("no auth"):
            headers = {"host": "localhost", "Content-Type": "application/json"}
            response = self.app.get(path, headers)
            response.raise_for_status()
            result = json.loads(response.body)
            collections = result.get("collections")
            self.assertIsNotNone(collections)
            ids = [collection.get("id") for collection in collections]
            self.assertIn(public_owned, ids)
            self.assertIn(public_not_owned, ids)
            self.assertNotIn(private_owned, ids)
            self.assertNotIn(private_not_owned, ids)

        with self.subTest("auth"):
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.get(path, headers)
            response.raise_for_status()
            result = json.loads(response.body)
            collections = result.get("collections")
            self.assertIsNotNone(collections)
            ids = [collection.get("id") for collection in collections]
            self.assertIn(public_owned, ids)
            self.assertIn(public_not_owned, ids)
            self.assertIn(private_owned, ids)
            self.assertNotIn(private_not_owned, ids)
