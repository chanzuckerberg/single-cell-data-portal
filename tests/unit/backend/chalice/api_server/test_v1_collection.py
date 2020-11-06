import itertools
import json
import sys
import os
import unittest
from datetime import datetime
from multiprocessing import Process

from furl import furl

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
            self.assertListEqual(sorted(collection.keys()), ["created_at", "id"])
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
            self.assertListEqual(sorted(link.keys()), ["name", "type", "url"])

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
                }
            ],
            "description": "test_description",
            "id": "test_collection_id",
            "links": [
                {"type": "RAW_DATA", "name": "test_link_name", "url": "test_url"},
                {"type": "OTHER", "name": "test_summary_link_name", "url": "test_summary_url"},
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

