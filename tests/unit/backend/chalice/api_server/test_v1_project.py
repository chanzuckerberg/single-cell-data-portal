import json
import unittest
from datetime import datetime

from furl import furl

from backend.corpora.common.corpora_orm import ProjectStatus
from backend.corpora.common.entities import Project, Dataset
from tests.unit.backend.chalice.api_server import BaseAPITest
from tests.unit.backend.utils import BogusProjectParams, BogusDatasetParams


class TestProject(BaseAPITest, unittest.TestCase):
    def test__list_project__ok(self):
        path = "/v1/project"
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
            self.assertIn(expected_id, [p["id"] for p in actual_body["projects"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(None, actual_body.get("from_date"))

        with self.subTest("from_date"):
            test_url = furl(path=path, query_params={"from_date": from_date})
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["projects"]])
            self.assertEqual(None, actual_body.get("to_date"))
            self.assertEqual(actual_body["from_date"], from_date)

        with self.subTest("to_date"):
            test_url = furl(path=path, query_params={"to_date": to_date})
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
            self.assertIn(expected_id, [p["id"] for p in actual_body["projects"]])
            self.assertEqual(to_date, actual_body["to_date"])
            self.assertEqual(None, actual_body.get("from_date"))

        with self.subTest("from_date->to_date"):
            test_url = furl(path=path, query_params={"from_date": from_date, "to_date": to_date})
            response = self.app.get(test_url.url, headers=headers)
            response.raise_for_status()
            actual_body = json.loads(response.body)
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
                    "development_stage": "test_development_stage",
                    "development_stage_ontology": "test_development_stage_ontology",
                    "disease": "test_disease",
                    "disease_ontology": "test_disease_ontology",
                    "ethnicity": "test_ethnicity",
                    "ethnicity_ontology": "test_ethnicity_ontology",
                    "id": "test_dataset_id",
                    "name": "test_dataset_name",
                    "organism": "test_organism",
                    "organism_ontology": "test_organism_ontology",
                    "preprint_doi": {"title": "test_preprint_doi"},
                    "project_id": "test_project_id",
                    "project_status": "LIVE",
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
            "owner": {"email": "test_email", "id": "test_user_id", "name": "test_user",},  # noqa
            "processing_state": "NA",
            "s3_bucket_key": "test_s3_bucket",
            "status": "LIVE",
            "validation_state": "NOT_VALIDATED",
        }

        test_url = furl(path="/v1/project/test_project_id")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        response.raise_for_status()
        actual_body = self.remove_timestamps(json.loads(response.body))
        actual_json_body = json.dumps(actual_body, sort_keys=True)
        expected_json_body = json.dumps(expected_body)
        self.assertEqual(actual_json_body, expected_json_body)

    def test__get_project_uuid__403_not_found(self):
        """Verify the test project exists and the expected fields exist."""
        test_url = furl(path="/v1/project/AAAA-BBBB-CCCC-DDDD")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, response.status_code)
        self.assertIn("X-AWS-REQUEST-ID", response.headers.keys())

    def test__delete_project_uuid__ok(self):
        expected_name = "test__delete_project_uuid__ok"
        test_project = Project.create(**BogusProjectParams.get(name=expected_name, status=ProjectStatus.LIVE.name))

        # check if it exists
        test_url = furl(path=f"/v1/project/{test_project.id}")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        response.raise_for_status()

        # delete
        test_url = furl(path=f"/v1/project/{test_project.id}")
        response = self.app.delete(test_url.url, headers=dict(host="localhost"))
        response.raise_for_status()

        # check if deleted
        test_url = furl(path=f"/v1/project/{test_project.id}")
        response = self.app.get(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, response.status_code)

    def test__delete_project_uuid__403_not_found(self):
        """Verify the test project exists and the expected fields exist."""
        test_url = furl(path="/v1/project/AAAA-BBBB-CCCC-DDDD")
        response = self.app.delete(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, response.status_code)


class TestProjectDataset(BaseAPITest, unittest.TestCase):
    def test__delete_dataset__ok(self):
        test_project = Project.create(**BogusProjectParams.get(status=ProjectStatus.EDIT.name))
        test_dataset = Dataset.create(
            **BogusDatasetParams.get(project_id=test_project.id, project_status=test_project.status)
        )
        expected_dataset_id = test_dataset.id

        # Verify Dataset Exists
        submission_url = furl(path=f"/v1/submission/{test_project.id}")
        resp = self.app.get(submission_url.url, headers=dict(host="localhost"))
        resp_json = json.loads(resp.body)

        actual_dataset_ids = [dataset["id"] for dataset in resp_json["datasets"]]
        self.assertIn(expected_dataset_id, actual_dataset_ids)

        # Delete the dataset
        dataset_url = furl(path=f"/v1/submission/{test_project.id}/dataset/{expected_dataset_id}")
        resp = self.app.delete(dataset_url.url, headers=dict(host="localhost"))
        self.assertEqual(202, resp.status_code)

        # Verify Dataset Deleted
        resp = self.app.get(submission_url.url, headers=dict(host="localhost"))
        resp_json = json.loads(resp.body)
        expected_dataset_id = test_dataset.id
        actual_dataset_ids = [dataset["id"] for dataset in resp_json["datasets"]]
        self.assertNotIn(expected_dataset_id, actual_dataset_ids)

    def test__delete_dataset__403_not_found(self):
        # Delete the dataset
        dataset_url = furl(path=f"/v1/submission/test_project_id/dataset/AAAA-BBBB-CCCC-DDDD")
        resp = self.app.delete(dataset_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, resp.status_code)

    def test__delete_dataset__403_not_in_project(self):
        test_project = Project.create(**BogusProjectParams.get(status=ProjectStatus.EDIT.name))
        test_dataset = Dataset.create(
            **BogusDatasetParams.get(project_id=test_project.id, project_status=test_project.status)
        )
        expected_dataset_id = test_dataset.id

        # Verify Dataset Exists
        submission_url = furl(path=f"/v1/submission/{test_project.id}")
        resp = self.app.get(submission_url.url, headers=dict(host="localhost"))
        resp_json = json.loads(resp.body)

        actual_dataset_ids = [dataset["id"] for dataset in resp_json["datasets"]]
        self.assertIn(expected_dataset_id, actual_dataset_ids)

        # Delete the dataset
        dataset_url = furl(path=f"/v1/submission/test_project_id/dataset/{expected_dataset_id}")
        resp = self.app.delete(dataset_url.url, headers=dict(host="localhost"))
        self.assertEqual(403, resp.status_code)

        # Verify Dataset Exists
        submission_url = furl(path=f"/v1/submission/{test_project.id}")
        resp = self.app.get(submission_url.url, headers=dict(host="localhost"))
        resp_json = json.loads(resp.body)

        actual_dataset_ids = [dataset["id"] for dataset in resp_json["datasets"]]
        self.assertIn(expected_dataset_id, actual_dataset_ids)
