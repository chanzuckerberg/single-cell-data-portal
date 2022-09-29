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
from backend.corpora.common.utils.corpora_constants import CorporaConstants
from backend.corpora.common.utils.db_helpers import processing_status_updater
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest, get_cxguser_token
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
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
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
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(403, response.status_code)

    def test__minimal_status__ok(self):
        dataset = self.generate_dataset(
            self.session, processing_status={"upload_status": "WAITING", "upload_progress": 0.0}
        )
        processing_status_id = dataset.processing_status.id
        test_url = furl(path=f"/dp/v1/datasets/{dataset.id}/status")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
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
            self.assertEqual(json.loads(response.data)["upload_status"], status.name)

    def test__get_all_datasets_for_index(self):
        test_dataset_id = "test_dataset_id_for_index"
        dataset = self.generate_dataset(
            self.session,
            id=test_dataset_id,
            cell_count=42,
            mean_genes_per_cell=0.05,
            published_at=datetime.now(),
            revised_at=datetime.now(),
        )
        self.generate_dataset(self.session, id="test_dataset_id_for_index_tombstone", tombstone=True)
        self.generate_dataset(
            self.session,
            id="test_dataset_id_for_index_private",
            collection_id="test_collection_id_revision",
        )
        test_url = furl(path="/dp/v1/datasets/index")
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        ids = [d["id"] for d in body]
        self.assertIn("test_dataset_id_for_index", ids)
        self.assertNotIn("test_dataset_id_for_index_tombstone", ids)
        self.assertNotIn("test_dataset_id_for_index_private", ids)

        actual_dataset = None
        for d in body:
            if d["id"] == test_dataset_id:
                actual_dataset = d
        self.assertIsNotNone(actual_dataset)

        self.assertEqual(actual_dataset["id"], dataset.id)
        self.assertEqual(actual_dataset["name"], dataset.name)
        self.assertNotIn("description", actual_dataset)
        self.assertEqual(actual_dataset["collection_id"], dataset.collection_id)
        self.assertEqual(actual_dataset["assay"], dataset.assay)
        self.assertEqual(actual_dataset["tissue"], dataset.tissue)
        self.assertEqual(actual_dataset["disease"], dataset.disease)
        self.assertEqual(actual_dataset["sex"], dataset.sex)
        self.assertEqual(actual_dataset["ethnicity"], dataset.ethnicity)
        self.assertEqual(actual_dataset["organism"], dataset.organism)
        self.assertEqual(actual_dataset["development_stage"], dataset.development_stage)
        self.assertEqual(actual_dataset["cell_count"], dataset.cell_count)
        self.assertEqual(actual_dataset["cell_type"], dataset.cell_type)
        self.assertEqual(actual_dataset["is_primary_data"], dataset.is_primary_data.name)
        self.assertEqual(actual_dataset["mean_genes_per_cell"], dataset.mean_genes_per_cell)
        self.assertEqual(actual_dataset["explorer_url"], dataset.explorer_url)
        self.assertEqual(actual_dataset["published_at"], dataset.published_at.timestamp())
        self.assertEqual(actual_dataset["revised_at"], dataset.revised_at.timestamp())

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

    def test__enrich_tissue_with_ancestors_expands_correctly(self):
        dataset = {"tissue": [{"ontology_term_id": "UBERON:0002048", "label": "Test"}]}
        Dataset.enrich_tissue_with_ancestors(dataset)
        self.assertIn("tissue_ancestors", dataset)
        self.assertEqual(
            dataset["tissue_ancestors"],
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

    def test__enrich_tissue_with_ancestors_empty_key_ok(self):
        dataset = {}
        Dataset.enrich_tissue_with_ancestors(dataset)
        self.assertEqual(dataset, {})

    def test__enrich_tissue_with_ancestors_missing_key_ok(self):
        dataset = {"tissue": [{"ontology_term_id": "UBERON:non_existent", "label": "Test"}]}
        Dataset.enrich_tissue_with_ancestors(dataset)
        self.assertNotIn("tissue_ancestors", dataset)

    def test__enrich_cell_type_with_ancestors_expands_correctly(self):
        dataset = {"cell_type": [{"ontology_term_id": "CL:0000738", "label": "Test"}]}
        Dataset.enrich_cell_type_with_ancestors(dataset)
        self.assertIn("cell_type_ancestors", dataset)
        self.assertEqual(
            dataset["cell_type_ancestors"],
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

    def test__enrich_cell_type_with_ancestors_empty_key_ok(self):
        dataset = {}
        Dataset.enrich_cell_type_with_ancestors(dataset)
        self.assertEqual(dataset, {})

    def test__enrich_cell_type_with_ancestors_missing_key_ok(self):
        dataset = {"cell_type": [{"ontology_term_id": "CL:non_existent", "label": "Test"}]}
        Dataset.enrich_cell_type_with_ancestors(dataset)
        self.assertNotIn("cell_type_ancestors", dataset)

    def test__get_all_datasets_for_index_with_ontology_expansion(self):
        test_dataset_id = "test_dataset_id_for_index_2"
        dataset = self.generate_dataset(
            self.session,
            id=test_dataset_id,
            cell_count=42,
            development_stage=[{"ontology_term_id": "HsapDv:0000008", "label": "Test"}],
            tissue=[{"ontology_term_id": "UBERON:0002048", "label": "Test"}],
            cell_type=[{"ontology_term_id": "CL:0000738", "label": "Test"}],
            published_at=datetime.now(),
            revised_at=datetime.now(),
        )

        test_url = furl(path="/dp/v1/datasets/index")

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.get(test_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        actual_dataset = None
        for d in body:
            if d["id"] == test_dataset_id:
                actual_dataset = d
        self.assertIsNotNone(actual_dataset)

        self.assertEqual(actual_dataset["development_stage"], dataset.development_stage)
        self.assertEqual(
            actual_dataset["development_stage_ancestors"],
            ["HsapDv:0000008", "HsapDv:0000006", "HsapDv:0000002", "HsapDv:0000045", "HsapDv:0000001"],
        )

        self.assertEqual(actual_dataset["tissue"], dataset.tissue)
        self.assertEqual(
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

        self.assertEqual(actual_dataset["cell_type"], dataset.cell_type)
        self.assertEqual(
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
        artifact_2 = dict(
            filename=CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME,
            filetype=DatasetArtifactFileType.H5AD,
            user_submitted=True,
            s3_uri="s3://mock-bucket/raw.h5ad",
        )
        dataset = self.generate_dataset(self.session, id="test_dataset", artifacts=[artifact_0, artifact_1, artifact_2])

        test_url = furl(path=f"/dp/v1/datasets/{dataset.id}/assets")
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
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # Test while uploading
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

    def test__cancel_dataset_download__dataset_does_not_exist(self):
        test_url = "/dp/v1/datasets/missing_dataset_id"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__delete_uploaded_dataset__ok(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}

        # check dataset in collection
        collection_url = furl(path=f"/dp/v1/collections/{collection.id}")
        collection_url.add(query_params=dict(visibility=CollectionVisibility.PRIVATE.name))
        response = self.app.get(collection_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)
        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertIn(dataset.id, dataset_ids)

        # delete dataset
        test_url = f"/dp/v1/datasets/{dataset.id}"
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)

        # check dataset no longer returned in collection
        response = self.app.get(collection_url.url, headers=headers)
        self.assertEqual(200, response.status_code)
        body = json.loads(response.data)

        dataset_ids = [dataset["id"] for dataset in body["datasets"]]
        self.assertNotIn(dataset.id, dataset_ids)

    def test__call_delete_dataset__twice(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"

        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)

        self.assertEqual(response.status_code, 202)

        # delete again
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__get_deleted_dataset_status__returns_403(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id"
        )
        processing_status = {"upload_status": UploadStatus.UPLOADED, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}/status"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(200, response.status_code)
        actual_body = json.loads(response.data)
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
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(405, response.status_code)
        self.assertEqual("Cannot delete a public Dataset", json.loads(response.data)["detail"])

    def test__cancel_dataset_download__user_not_collection_owner(self):
        collection = self.generate_collection(self.session, owner="someone_else")
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 403)

    def test__cancel_dataset_download__user_not_logged_in(self):
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 401)

    def test__dataset_meta__ok(self):
        public_collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="test_user_id"
        )
        private_collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="someone_else"
        )
        headers = {"host": "localhost", "Content-Type": "application/json"}

        with self.subTest("dataset is public"):
            test_uri_0 = "some_uri_0"
            artifact_params_0 = dict(
                filename="filename_x",
                filetype=DatasetArtifactFileType.CXG,
                user_submitted=True,
                s3_uri=test_uri_0,
            )
            public_dataset = self.generate_dataset(
                self.session, collection=public_collection, explorer_url="test_url_0", artifacts=[artifact_params_0]
            )
            test_url_public = f"/dp/v1/datasets/meta?url={public_dataset.explorer_url}"

            response = self.app.get(test_url_public, headers)
            self.assertEqual(response.status_code, 200)

            expected_identifiers = {
                "s3_uri": test_uri_0,
                "dataset_id": public_dataset.id,
                "collection_id": public_collection.id,
                "collection_visibility": public_collection.visibility.name,
                "tombstoned": False,
            }

            self.assertEqual(json.loads(response.data), expected_identifiers)

        with self.subTest("dataset is private"):
            test_uri_1 = "some_uri_1"
            artifact_params_1 = dict(
                filename="filename_x",
                filetype=DatasetArtifactFileType.CXG,
                user_submitted=True,
                s3_uri=test_uri_1,
            )
            private_dataset = self.generate_dataset(
                self.session, collection=private_collection, explorer_url="test_url_1", artifacts=[artifact_params_1]
            )
            test_url_private = f"/dp/v1/datasets/meta?url={private_dataset.explorer_url}"
            expected_identifiers = {
                "s3_uri": test_uri_1,
                "dataset_id": private_dataset.id,
                "collection_id": private_collection.id,
                "collection_visibility": private_collection.visibility.name,
                "tombstoned": False,
            }

            response = self.app.get(test_url_private, headers)
            self.assertEqual(response.status_code, 200)

            self.assertEqual(json.loads(response.data), expected_identifiers)

        with self.subTest("dataset does not have an associated cxg artifact"):
            dataset_without_artifacts = self.generate_dataset(
                self.session, collection=public_collection, explorer_url="test_url_2"
            )
            test_url_no_cxg_artifact = f"/dp/v1/datasets/meta?url={dataset_without_artifacts.explorer_url}"
            expected_identifiers = {
                "s3_uri": None,
                "dataset_id": dataset_without_artifacts.id,
                "collection_id": public_collection.id,
                "collection_visibility": "PUBLIC",
                "tombstoned": False,
            }

            response = self.app.get(test_url_no_cxg_artifact, headers)
            self.assertEqual(response.status_code, 200)

            self.assertEqual(json.loads(response.data), expected_identifiers)

    def test__dataset_meta__404(self):
        headers = {"host": "localhost", "Content-Type": "application/json"}
        test_url_404 = "/dp/v1/datasets/meta?url=not_real"

        response = self.app.get(test_url_404, headers)
        self.assertEqual(response.status_code, 404)


class TestDatasetGenesetLinkageUpdates(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()

    def tearDown(self):
        super().tearDown()

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
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
        response = self.app.post(test_url, headers=headers, data=json.dumps(data))
        self.assertEqual(len(json.loads(response.data)), 4)
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
            test_url = f"/dp/v1/datasets/{generate_id()}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 403)

        with self.subTest("collection public"):
            test_url = f"/dp/v1/datasets/{dataset_0.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 403)

        with self.subTest("user not collection owner"):
            test_url = f"/dp/v1/datasets/{dataset_1.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
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
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 404)

        with self.subTest("remove list references genesets that do not belong to the collection"):
            data = {"add": [], "remove": [geneset1.id]}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 404)

        with self.subTest("add list references genesets already linked to the dataset"):
            data = {"add": [geneset0.id], "remove": []}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 404)

        with self.subTest("remove list references genesets that are not currently linked to the dataset"):
            data = {"add": [], "remove": [geneset2.id]}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 404)

        with self.subTest("request body missing remove list"):
            data = {"remove": [geneset0.id]}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 400)

        with self.subTest("request body missing add list"):
            data = {"add": [geneset2.id]}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 400)

        with self.subTest("no data"):
            data = {}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 400)

        with self.subTest("empty data"):
            data = {"add": [], "remove": []}
            test_url = f"/dp/v1/datasets/{dataset.id}/gene_sets"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token()}
            response = self.app.post(test_url, headers=headers, data=json.dumps(data))
            self.assertEqual(response.status_code, 202)


class TestDatasetCurators(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()

    def tearDown(self):
        super().tearDown()

    def test__get_status__200_for_non_owned_dataset_as_super_curator(self):
        collection = self.generate_collection(self.session, owner="someone_else")
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}/status"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token("super")}
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(200, response.status_code)

    def test__cancel_dataset_download__202_user_not_collection_owner_as_super_curator(self):
        collection = self.generate_collection(self.session, owner="someone_else")
        processing_status = {"upload_status": UploadStatus.WAITING, "upload_progress": 0.0}
        dataset = self.generate_dataset(self.session, collection=collection, processing_status=processing_status)
        test_url = f"/dp/v1/datasets/{dataset.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token("super")}
        response = self.app.delete(test_url, headers=headers)
        self.assertEqual(response.status_code, 202)
