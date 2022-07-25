import json
import os
import unittest

from furl import furl

from datetime import datetime

from backend.corpora.common.corpora_orm import (
    UploadStatus,
    CollectionVisibility,
    generate_id,
    DatasetArtifactFileType,
)
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.utils.db_helpers import processing_status_updater
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest, BasicAuthAPITestCurator
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestDataset(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()

    def tearDown(self):
        super().tearDown()

    def test__post_dataset_asset__OK(self):
        bucket = self.CORPORA_TEST_CONFIG["bucket_name"]
        s3_file_name = "test_s3_uri.h5ad"
        content = "Hello world!"
        self.create_s3_object(s3_file_name, bucket, content=content)

        expected_body = dict(dataset_id="test_dataset_id", file_name="test_filename", file_size=len(content))
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        presign_url = actual_body.pop("presigned_url")
        self.assertIsNotNone(presign_url)
        self.assertEqual(expected_body, actual_body)

    @unittest.skipIf(os.getenv("IS_DOCKER_DEV"), "This is currently TODO in docker")
    def test__post_dataset_asset__file_SERVER_ERROR(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(500, response.status_code)
        body = json.loads(response.data)
        self.assertEqual("An internal server error has occurred. Please try again later.", body["detail"])

    def test__post_dataset_asset__dataset_NOT_FOUND(self):
        test_url = furl(path="/dp/v1/datasets/test_user_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.data)
        self.assertEqual("'dataset/test_user_id' not found.", body["detail"])
        print(body)

    def test__post_dataset_asset__asset_NOT_FOUND(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/fake_asset")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.data)
        self.assertEqual("'dataset/test_dataset_id/asset/fake_asset' not found.", body["detail"])

    def test__get_status__ok(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        expected_body = {
            "cxg_status": "NA",
            "rds_status": "NA",
            "h5ad_status": "NA",
            "processing_status": "PENDING",
            "dataset_id": "test_dataset_id",
            "id": "test_dataset_processing_status_id",
            "upload_progress": 0.4444444444444444,
            "upload_status": "UPLOADING",
            "validation_status": "NA",
        }
        self.assertEqual(expected_body, actual_body)

    def test__get_status__403(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id_not_owner/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(403, response.status_code)

    def test__minimal_status__ok(self):
        id, _ = self.generate_dataset(processing_status={"upload_status": "WAITING", "upload_progress": 0.0})
        test_url = furl(path=f"/dp/v1/datasets/{id}/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        expected_body = {
            "upload_progress": 0.0,
            "upload_status": "WAITING",
        }
        self.assertEqual(expected_body, actual_body)

        for status in UploadStatus:
            processing_status = {"upload_status": status, "upload_progress": 0.0}
            self.update_processing(id, processing_status)
            response = self.app.get(test_url.url, headers=headers)
            self.assertEqual(json.loads(response.data)["upload_status"], status.name)

    def test__get_all_datasets_for_index(self):
        coll1, _ = self.generate_collection()
        self.publish_collection(coll1)
        coll2, _ = self.generate_collection()
        id, dataset = self.generate_dataset(coll1, cell_count=42)
        id2, _ = self.generate_dataset(coll2, cell_count=43)
        test_url = furl(path="/dp/v1/datasets/index")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        counts = [d["cell_count"] for d in body]
        self.assertIn(42, counts)
        self.assertNotIn(43, counts)

        actual_dataset = None
        for d in body:
            if d["cell_count"] == 42:
                actual_dataset = d
        self.assertIsNotNone(actual_dataset)

        self.assertDictEqual(self.remove_timestamps(actual_dataset), self.remove_timestamps(dataset))

    def test__enrich_development_stage_with_ancestors_expands_correctly(self):
        dataset = {"development_stage": [{"ontology_term_id": "HsapDv:0000008", "label": "Test"}]}
        Dataset.enrich_development_stage_with_ancestors(dataset)
        self.assertIn("development_stage_ancestors", dataset)
        self.assertEqual(
            dataset["development_stage_ancestors"],
            ["HsapDv:0000008", "HsapDv:0000006", "HsapDv:0000002", "HsapDv:0000045", "HsapDv:0000001"],
        )

    def test__enrich_development_stage_with_ancestors_empty_key_ok(self):
        dataset = {}
        Dataset.enrich_development_stage_with_ancestors(dataset)
        self.assertEqual(dataset, {})

    def test__enrich_development_stage_with_ancestors_missing_key_ok(self):
        dataset = {"development_stage": [{"ontology_term_id": "HsapDv:non_existant", "label": "Test"}]}
        Dataset.enrich_development_stage_with_ancestors(dataset)
        self.assertNotIn("development_stage_ancestors", dataset)

    def test__get_all_datasets_for_index_with_ontology_expansion(self):
        dataset, _ = self.generate_dataset(
            self.session,
            cell_count=42,
            development_stage=[{"ontology_term_id": "HsapDv:0000008", "label": "Test"}],
        )

        test_url = furl(path="/dp/v1/datasets/index")

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        actual_dataset = None
        for d in body:
            if d["id"] == dataset:
                actual_dataset = d
        self.assertIsNotNone(actual_dataset)

        self.assertEqual(actual_dataset["development_stage"], dataset['development_stage'])
        self.assertEqual(
            actual_dataset["development_stage_ancestors"],
            ["HsapDv:0000008", "HsapDv:0000006", "HsapDv:0000002", "HsapDv:0000045", "HsapDv:0000001"],
        )

    def test__get_dataset_assets(self):
        artifact_0 = dict(
            filename="filename_0",
            filetype=DatasetArtifactFileType.CXG,
            user_submitted=True,
            s3_uri="s3://mock-bucket/mock-key.cxg",
        )
        artifact_1 = dict(
            filename="filename_1",
            filetype=DatasetArtifactFileType.H5AD,
            user_submitted=True,
            s3_uri="s3://mock-bucket/mock-key.h5ad",
        )
        coll_id, _ = self.generate_collection()
        id, dataset = self.generate_dataset(coll_id, artifacts=[artifact_0, artifact_1])

        test_url = furl(path=f"/dp/v1/datasets/{id}/assets")
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        self.assertIn("assets", body)
        assets = body["assets"]
        self.assertEqual(len(assets), 2)
        self.assertEqual(assets[0]["s3_uri"], "s3://mock-bucket/mock-key.cxg")
        self.assertEqual(assets[1]["s3_uri"], "s3://mock-bucket/mock-key.h5ad")

    def test__cancel_dataset_download__ok(self):
        # Test pre upload
        collection, _ = self.generate_collection()
        processing_status = {"upload_status": UploadStatus.WAITING.name, "upload_progress": 0.0}
        dataset, _ = self.generate_dataset(collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # Test while uploading
        processing_status = {"upload_status": UploadStatus.UPLOADING.name, "upload_progress": 10.0}
        dataset, _ = self.generate_dataset(collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

    def test__cancel_dataset_download__dataset_does_not_exist(self):
        test_url = "/dp/v1/datasets/missing_dataset_id"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__delete_uploaded_dataset__ok(self):
        coll_id, _ = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.UPLOADED.name, "upload_progress": 0.0}
        dataset_id, _ = self.generate_dataset(coll_id, processing_status=processing_status)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        # check dataset in collection
        collection_url = furl(path=f"/dp/v1/collections/{coll_id}")
        response = self.app.get(collection_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset_id, dataset_ids)

        # delete dataset
        test_url = f"/dp/v1/datasets/{dataset_id}?collection_id={id}"
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

        # check dataset no longer returned in collection
        response = self.app.get(collection_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertNotIn(dataset_id, dataset_ids)

    def test__call_delete_dataset__twice(self):
        id, _ = self.generate_collection()
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        dataset_id, _ = self.generate_dataset(id, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset_id}?collection_id={id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # delete again
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__get_deleted_dataset_status__returns_403(self):
        collection, _ = self.generate_collection(owner="test_user_id")
        processing_status = {"upload_status": UploadStatus.UPLOADED.name, "upload_progress": 0.0}
        dataset, _ = self.generate_dataset(collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset}/status"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
        expected_body = {
            "upload_progress": 0.0,
            "upload_status": "UPLOADED",
        }
        self.assertEqual(expected_body, actual_body)

        # delete the dataset
        self.app.delete(f"/dp/v1/datasets/{dataset}", headers=headers)

        response = self.app.get(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__delete_public_dataset_returns__405(self):
        collection, _ = self.generate_collection()
        self.publish_collection(collection)
        processing_status = {"upload_status": UploadStatus.UPLOADED.name, "upload_progress": 0.0}
        dataset, _ = self.generate_dataset(collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(405, response.status_code)
        self.assertEqual("Cannot delete a public Dataset", json.loads(response.data)["detail"])

    def test__cancel_dataset_download__user_not_collection_owner(self):
        collection, _ = self.generate_collection(owner="someone_else")
        processing_status = {"upload_status": UploadStatus.WAITING.name, "upload_progress": 0.0}
        dataset, _ = self.generate_dataset(collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__cancel_dataset_download__user_not_logged_in(self):
        processing_status = {"upload_status": UploadStatus.WAITING.name, "upload_progress": 0.0}
        collection, _ = self.generate_collection()
        dataset, _ = self.generate_dataset(collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset}"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 401)

    def test__dataset_meta__ok(self):
        public_id, public_collection = self.generate_collection(owner="test_user_id")
        self.publish_collection(public_id)
        private_id, private_collection = self.generate_collection(owner="someone_else")
        headers = {"host": "localhost", "Content-Type": "application/json"}

        with self.subTest("dataset is public"):
            test_uri_0 = "some_uri_0"
            artifact_params_0 = dict(
                filename="filename_x",
                filetype=DatasetArtifactFileType.CXG,
                user_submitted=True,
                s3_uri=test_uri_0,
            )
            public_dataset_id, public_dataset = self.generate_dataset(
                public_id, explorer_url="test_url_0", artifacts=[artifact_params_0]
            )
            test_url_public = f"/dp/v1/datasets/meta?url={public_dataset['explorer_url']}"

            response = self.app.get(test_url_public, headers)
            self.assertEqual(response.status_code, 200)

            expected_identifiers = {
                "s3_uri": test_uri_0,
                "dataset_id": public_dataset_id,
                "collection_id": public_id,
                "collection_visibility": public_collection['visibility'],
            }

            self.assertEqual(json.loads(response.data), expected_identifiers)

        with self.subTest("dataset is private"):
            test_uri_1 = "some_uri_1"
            artifact_params_1 = dict(
                filename="filename_x",
                filetype=DatasetArtifactFileType.CXG.name,
                user_submitted=True,
                s3_uri=test_uri_1,
            )
            private_dataset_id, private_dataset = self.generate_dataset(
                private_id, explorer_url="test_url_1", artifacts=[artifact_params_1]
            )
            test_url_private = f"/dp/v1/datasets/meta?url={private_dataset['explorer_url']}"
            expected_identifiers = {
                "s3_uri": test_uri_1,
                "dataset_id": private_dataset_id,
                "collection_id": private_id,
                "collection_visibility": private_collection['visibility'],
            }

            response = self.app.get(test_url_private, headers)
            self.assertEqual(response.status_code, 200)

            self.assertEqual(json.loads(response.data), expected_identifiers)

        with self.subTest("dataset does not have an associated cxg artifact"):
            dataset_without_artifacts_id, dataset_without_artifacts = self.generate_dataset(
                public_id, explorer_url="test_url_2"
            )
            test_url_no_cxg_artifact = f"/dp/v1/datasets/meta?url={dataset_without_artifacts['explorer_url']}"
            expected_identifiers = {
                "s3_uri": None,
                "dataset_id": dataset_without_artifacts_id,
                "collection_id": public_id,
                "collection_visibility": "PUBLIC",
            }

            response = self.app.get(test_url_no_cxg_artifact, headers)
            self.assertEqual(response.status_code, 200)

            self.assertEqual(json.loads(response.data), expected_identifiers)

    def test__dataset_meta__404(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}
        test_url_404 = "/dp/v1/datasets/meta?url=not_real"

        response = self.app.get(test_url_404, headers)
        self.assertEqual(response.status_code, 404)


class TestDatasetCurators(BasicAuthAPITestCurator, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()

    def tearDown(self):
        super().tearDown()

    def test__get_status__200_for_non_owned_dataset(self):
        id, _ = self.generate_collection(owner="someone_else")
        processing_status = {"upload_status": UploadStatus.WAITING.name, "upload_progress": 0.0}
        dataset_id, _ = self.generate_dataset(id, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset_id}/status?collection_id={id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(200, response.status_code)

    def test__cancel_dataset_download__202_user_not_collection_owner(self):
        id, _ = self.generate_collection(owner="someone_else")
        processing_status = {"upload_status": UploadStatus.WAITING.name, "upload_progress": 0.0}
        dataset_id, _ = self.generate_dataset(id, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset_id}?collection_id={id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)
