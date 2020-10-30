import json
import unittest
from datetime import datetime

from furl import furl

from backend.corpora.common.corpora_orm import ProjectStatus
from backend.corpora.common.entities import Project
from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.utils import BogusProjectParams


class TestProject(BaseAPITest, unittest.TestCase):
    def validate_projects_response_structure(self, body):
        self.assertIn("projects", body)
        self.assertTrue(all(k in ["projects", "from_date", "to_date"] for k in body))

        for project in body["projects"]:
            self.assertListEqual(sorted(project.keys()), ["created_at", "id"])
            self.assertGreaterEqual(datetime.fromtimestamp(project["created_at"]).year, 1969)

    def validate_project_uuid_response_structure(self, body):
        required_keys = [
            "name",
            "description",
            "id",
            "s3_bucket_key",
            "status",
            "processing_state",
            "validation_state",
            "links",
            "attestation",
            "datasets",
            "created_at",
            "updated_at",
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
                "source_data_location",
                "revision",
                "dataset_deployments",
                "dataset_assets",
                "preprint_doi",
                "publication_doi",
                "created_at",
                "updated_at",
                "project_id",
                "project_status",
            ]
            self.assertListEqual(sorted(dataset.keys()), sorted(required_keys))

    def test__list_project__ok(self):
        path = "/dp/v1/project"
        headers = dict(host="localhost")

        from_date = int(datetime.fromtimestamp(60).timestamp())
        creation_time = 70
        to_date = int(datetime.fromtimestamp(80).timestamp())

        test_project = Project.create(
            **BogusProjectParams.get(status=ProjectStatus.LIVE.name, created_at=datetime.fromtimestamp(creation_time)),
        )
        expected_id = test_project.id

        with self.subTest("No Parameters"):
            test_url = furl(path=path)
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.validate_projects_response_structure(actual_body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["projects"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(None, actual_body.get("from_date"))

        with self.subTest("from_date"):
            test_url = furl(path=path, query_params={"from_date": from_date})
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.validate_projects_response_structure(actual_body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["projects"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(actual_body["from_date"], from_date)

        with self.subTest("to_date"):
            test_url = furl(path=path, query_params={"to_date": to_date})
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.validate_projects_response_structure(actual_body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["projects"]])
            self.assertEqual(to_date, actual_body["to_date"])
            self.assertEqual(None, actual_body.get("from_date"))

        with self.subTest("from_date->to_date"):
            test_url = furl(path=path, query_params={"from_date": from_date, "to_date": to_date})
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.validate_projects_response_structure(actual_body)
            self.assertEqual(expected_id, actual_body["projects"][0]["id"])
            self.assertEqual(creation_time, actual_body["projects"][0]["created_at"])
            self.assertEqual(from_date, actual_body["from_date"])
            self.assertEqual(to_date, actual_body["to_date"])

    def test__get_project_uuid__ok(self):
        """Verify the test project exists and the expected fields exist."""
        expected_body = {
            "attestation": {"needed": False, "tc_uri": "test_tc_uri"},
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
                            "environment": "test",
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
                    "preprint_doi": "test_preprint_doi",
                    "project_id": "test_project_id",
                    "project_status": "LIVE",
                    "publication_doi": "test_publication_doi",
                    "revision": 0,
                    "sex": ["test_sex", "test_sex2"],
                    "tissue": [{"label": "test_tissue", "ontology_term_id": "test_obo"}],
                    "source_data_location": "test_source_data_location",
                }
            ],
            "description": "test_description",
            "id": "test_project_id",
            "links": [
                {"type": "RAW_DATA", "name": "test_link_name", "url": "test_url"},
                {"type": "SUMMARY", "name": "test_summary_link_name", "url": "test_summary_url"},
            ],
            "name": "test_project",
            "processing_state": "NA",
            "s3_bucket_key": "test_s3_bucket",
            "status": "LIVE",
            "validation_state": "NOT_VALIDATED",
        }

        test_url = furl(path="/dp/v1/project/test_project_id")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        response.raise_for_status()
        self.validate_project_uuid_response_structure(json.loads(response.body))
        actual_body = self.remove_timestamps(json.loads(response.body))
        self.assertDictEqual(actual_body, expected_body)

    def test__get_project_uuid__403_not_found(self):
        """Verify the test project exists and the expected fields exist."""
        test_url = furl(path="/dp/v1/project/AAAA-BBBB-CCCC-DDDD")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, response.status_code)
        self.assertIn("X-AWS-REQUEST-ID", response.headers.keys())
