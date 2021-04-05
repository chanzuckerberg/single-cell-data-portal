import json
import os
import unittest
from requests import HTTPError

from furl import furl

from backend.corpora.common.corpora_orm import UploadStatus, CollectionVisibility, generate_uuid
from backend.corpora.common.utils.db_session import processing_status_updater
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestDataset(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        CorporaTestCaseUsingMockAWS.setUp(self)

    def tearDown(self):
        CorporaTestCaseUsingMockAWS.tearDown(self)

    def test__post_dataset_asset__OK(self):
        bucket = self.CORPORA_TEST_CONFIG["bucket_name"]
        s3_file_name = "test_s3_uri.h5ad"
        content = "Hello world!"
        self.create_s3_object(s3_file_name, bucket, content=content)

        expected_body = dict(dataset_id="test_dataset_id", file_name="test_filename", file_size=len(content))
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        response.raise_for_status()
        actual_body = json.loads(response.body)
        presign_url = actual_body.pop("presigned_url")
        self.assertIsNotNone(presign_url)
        self.assertEqual(expected_body, actual_body)

    @unittest.skipIf(os.getenv("IS_DOCKER_DEV"), "This is currently TODO in docker")
    def test__post_dataset_asset__file_SERVER_ERROR(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(500, response.status_code)
        body = json.loads(response.body)
        self.assertEqual("An internal server error has occurred. Please try again later.", body["detail"])

    def test__post_dataset_asset__dataset_NOT_FOUND(self):
        test_url = furl(path="/dp/v1/datasets/test_user_id/asset/test_dataset_artifact_id")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.body)
        self.assertEqual("'dataset/test_user_id' not found.", body["detail"])
        print(body)

    def test__post_dataset_asset__asset_NOT_FOUND(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/asset/fake_asset")
        response = self.app.post(test_url.url, headers=dict(host="localhost"))
        self.assertEqual(404, response.status_code)
        body = json.loads(response.body)
        self.assertEqual("'dataset/test_dataset_id/asset/fake_asset' not found.", body["detail"])

    def test__get_status__ok(self):
        test_url = furl(path="/dp/v1/datasets/test_dataset_id/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        response.raise_for_status()
        actual_body = json.loads(response.body)
        expected_body = {
            "conversion_anndata_status": "NA",
            "conversion_cxg_status": "NA",
            "conversion_loom_status": "NA",
            "conversion_rds_status": "NA",
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
        dataset = self.generate_dataset(
            self.session, processing_status={"upload_status": "WAITING", "upload_progress": 0.0}
        )
        processing_status_id = dataset.processing_status.id
        test_url = furl(path=f"/dp/v1/datasets/{dataset.id}/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url.url, headers=headers)
        response.raise_for_status()
        actual_body = json.loads(response.body)
        expected_body = {
            "dataset_id": dataset.id,
            "id": processing_status_id,
            "upload_progress": 0.0,
            "upload_status": "WAITING",
        }
        self.assertEqual(expected_body, actual_body)

        for status in UploadStatus:
            processing_status = {"upload_status": status, "upload_progress": 0.0}
            processing_status_updater(self.session, processing_status_id, processing_status)
            response = self.app.get(test_url.url, headers=headers)
            self.assertEqual(json.loads(response.body)["upload_status"], status.name)

    def test__cancel_dataset_download__ok(self):
        # Test pre upload
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # Test while uploading
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

    def test__cancel_dataset_download__dataset_does_not_exist(self):
        test_url = "/dp/v1/datasets/missing_dataset_id"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__delete_uploaded_dataset__ok(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        # check dataset in collection
        collection_url = furl(path=f"/dp/v1/collections/{collection.id}")
        collection_url.add(query_params=dict(visibility=CollectionVisibility.PRIVATE.name))
        response = self.app.get(collection_url.url, headers=headers)
        response.raise_for_status()
        body = json.loads(response.body)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset.id, dataset_ids)

        # delete dataset
        test_url = f"/dp/v1/datasets/{dataset.id}"
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

        # check dataset no longer returned in collection
        response = self.app.get(collection_url.url, headers=headers)
        response.raise_for_status()
        body = json.loads(response.body)

        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertNotIn(dataset.id, dataset_ids)

    def test__delete_dataset_can_safetly_be_called_twice(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # delete again
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

    def test__get_deleted_dataset_status__returns_403(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}/status"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.get(test_url, headers=headers)
        response.raise_for_status()
        actual_body = json.loads(response.body)
        expected_body = {
            "dataset_id": dataset.id,
            "id": dataset.processing_status.id,
            "upload_progress": 0.0,
            "upload_status": "UPLOADED",
        }
        self.assertEqual(expected_body, actual_body)

        # delete the dataset
        self.app.delete(f"/dp/v1/datasets/{dataset.id}", headers=headers)

        response = self.app.get(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__delete_public_dataset_returns__405(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name)
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 405)
        self.assertEqual(json.loads(response.body), "Can not delete a public dataset")

    def test__cancel_dataset_download__user_not_collection_owner(self):
        collection = self.generate_collection(self.session, owner="someone_else")
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__cancel_dataset_download__user_not_logged_in(self):
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 401)


class TestDatasetGenesetLinkageUpdates(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        CorporaTestCaseUsingMockAWS.setUp(self)

    def tearDown(self):
        CorporaTestCaseUsingMockAWS.tearDown(self)

    def test__dataset_gene_set_linkage_update__ok(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        dataset = self.generate_dataset(self.session, collection=collection)
        geneset0 = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        geneset1 = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        geneset2 = self.generate_geneset(self.session, collection=collection, dataset_ids=[dataset.id])
        geneset3 = self.generate_geneset(self.session, collection=collection)
        geneset4 = self.generate_geneset(self.session, collection=collection)

        links = dataset.genesets
        self.assertEqual(len(links), 3)
        link_ids = [x.id for x in links]
        self.assertIn(geneset0.id, link_ids)
        self.assertIn(geneset1.id, link_ids)
        self.assertIn(geneset2.id, link_ids)
        data = {"add": [geneset3.id, geneset4.id], "remove": [geneset2.id]}
        test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.post(test_url, headers, data=json.dumps(data))
        response.raise_for_status()
        self.assertEqual(len(json.loads(response.body)), 4)
        self.assertEqual(response.status_code, 202)

    def test__dataset_gene_set_linkage_update__403(self):
        collection_0 = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )
        collection_1 = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="someone_else"
        )
        dataset_0 = self.generate_dataset(self.session, collection=collection_0)
        dataset_1 = self.generate_dataset(self.session, collection=collection_1)

        with self.subTest("dataset does not exist"):
            data = {"add": [], "remove": []}
            test_url = f"/dp/v1/datasets/{generate_uuid()}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 403)

        with self.subTest("collection public"):
            test_url = f"/dp/v1/datasets/{dataset_0.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 403)

        with self.subTest("user not collection owner"):
            test_url = f"/dp/v1/datasets/{dataset_1.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 403)

    def test__dataset_gene_set_linkage_update__404(self):
        collection_0 = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        dataset = self.generate_dataset(self.session, collection=collection_0)
        collection_1 = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        geneset0 = self.generate_geneset(self.session, collection=collection_0, dataset_ids=[dataset.id])
        geneset1 = self.generate_geneset(self.session, collection=collection_1)
        geneset2 = self.generate_geneset(self.session, collection=collection_0)

        with self.subTest("add list references genesets that do not belong to the collection"):
            data = {"add": [geneset1.id], "remove": []}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 404)

        with self.subTest("remove list references genesets that do not belong to the collection"):
            data = {"add": [], "remove": [geneset1.id]}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 404)

        with self.subTest("add list references genesets already linked to the dataset"):
            data = {"add": [geneset0.id], "remove": []}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 404)

        with self.subTest("remove list references genesets that are not currently linked to the dataset"):
            data = {"add": [], "remove": [geneset2.id]}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 404)

        with self.subTest("request body missing remove list"):
            data = {"remove": [geneset0.id]}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 400)

        with self.subTest("request body missing add list"):
            data = {"add": [geneset2.id]}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 400)

        with self.subTest("no data"):
            data = {}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            with self.assertRaises(HTTPError):
                response.raise_for_status()
            self.assertEqual(response.status_code, 400)

        with self.subTest("empty data"):
            data = {"add": [], "remove": []}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            response = self.app.post(test_url, headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 202)
