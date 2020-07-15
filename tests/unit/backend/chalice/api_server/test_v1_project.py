import json
import unittest
from datetime import datetime

from furl import furl

from backend.corpora.common.corpora_orm import ProjectStatus
from backend.corpora.common.entities import Project
from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.utils import BogusProjectParams


class TestProject(BaseAPITest, unittest.TestCase):
    def test__list_project__ok(self):
        with self.subTest("No Parameters"):
            test_url = furl(path="/v1/project")
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            result_body = json.loads(response.body)
            self.assertIn("test_project_id", [p["id"] for p in result_body["projects"]])

        creation_time = 0

        test_project = Project.create(
            **BogusProjectParams.get(status=ProjectStatus.LIVE.name), created_at=datetime.fromtimestamp(creation_time)
        )
        test_id = test_project.id
        future_time = int(datetime.fromtimestamp(10).timestamp())
        with self.subTest("With to_date"):
            test_url = furl(path="/v1/project", query_params={"to_date": future_time})
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            result_body = json.loads(response.body)
            self.assertEqual(result_body["projects"][0]["id"], test_id)
            self.assertEqual(result_body["projects"][0]["created_at"], creation_time)
            self.assertEqual(result_body["to_date"], future_time)

        with self.subTest("With from_date"):
            test_url = furl(path="/v1/project", query_params={"from_date": future_time})
            response = self.app.get(test_url.url, headers=dict(host="localhost"))
            response.raise_for_status()
            result_body = json.loads(response.body)
            self.assertIn("test_project_id", [p["id"] for p in result_body["projects"]])
            self.assertEqual(result_body["from_date"], future_time)

    @staticmethod
    def remove_timestamps(json_response_body: dict) -> dict:
        def _remove_timestamps(jrb):
            jrb.pop("created_at", None)
            jrb.pop("updated_at", None)
            for value in jrb.values():
                if isinstance(value, dict):
                    _remove_timestamps(value)
                elif isinstance(value, list):
                    for list_value in value:
                        _remove_timestamps(list_value)
            return jrb

        return _remove_timestamps(json_response_body)

    def test__get_project_uuid__ok(self):
        """Verify the test project exists and the expected fields exist."""
        expected_response_body = {
            "attestation": {"needed": False, "tc_uri": "test_tc_uri"},
            "datasets": [
                {
                    "assay": "test_assay",
                    "assay_ontology": "test_assay_ontology",
                    "contributors": [
                        {
                            "email": "test_email",
                            "id": "test_contributor_id",
                            "institution": "test_institution",
                            "name": "test_contributor_name",
                        }
                    ],
                    "dataset_assets": [
                        {
                            "dataset_id": "test_dataset_id",
                            "filename": "test_filename",
                            "filetype": "H5AD",
                            "id": "test_dataset_artifact_id",
                            "s3_uri": "test_s3_uri",
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
                    "disease": "test_disease",
                    "disease_ontology": "test_disease_ontology",
                    "ethnicity": "test_ethnicity",
                    "ethnicity_ontology": "test_ethnicity_ontology",
                    "id": "test_dataset_id",
                    "name": "test_dataset_name",
                    "organism": "test_organism",
                    "organism_ontology": "test_organism_ontology",
                    "preprint_doi": {"title": "test_preprint_doi"},
                    "publication_doi": {"title": "test_publication_doi"},
                    "revision": 0,
                    "sex": "test_sex",
                    "source_data_location": "test_source_data_location",
                    "tissue": "test_tissue",
                    "tissue_ontology": "test_tissue_ontology",
                }
            ],
            "description": "test_description",
            "id": "test_project_id",
            "links": [{"type": "RAW_DATA", "url": "test_url"}],
            "name": "test_project",
            "owner": {"email": "test_email", "id": "test_user_id", "name": "test_user",},
            "processing_state": "NA",
            "s3_bucket_key": "test_s3_bucket",
            "status": "LIVE",
            "validation_state": "NOT_VALIDATED",
        }

        test_url = furl(path="/v1/project/test_project_id")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        response.raise_for_status()
        result_body = json.loads(response.body)
        result_body = self.remove_timestamps(result_body)
        result_json_body = json.dumps(result_body, sort_keys=True)
        expected_json_body = json.dumps(expected_response_body)
        self.assertEqual(result_json_body, expected_json_body)
