import unittest

from backend.common.corpora_orm import CollectionVisibility, IsPrimaryData, UploadStatus
from backend.layers.common.entities import DatasetStatusKey, DatasetUploadStatus
from tests.unit.backend.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.layers.common.base_api_test import DatasetStatusUpdate


class TestDeleteDataset(BaseAuthAPITest):
    def test__delete_dataset(self):
        auth_credentials = [
            (self.make_super_curator_header, "super", 202),
            (self.make_owner_header, "owner", 202),
            (None, "none", 401),
            (self.make_not_owner_header, "not_owner", 403),
        ]
        for auth, auth_description, expected_status_code in auth_credentials:
            with self.subTest(f"{auth_description} {expected_status_code}"):

                dataset = self.generate_dataset(
                    statuses=[DatasetStatusUpdate(DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)],
                    publish=False
                )

                test_url = f"/curation/v1/collections/{dataset.collection_version_id}/datasets/{dataset.dataset_version_id}"
                headers = auth() if callable(auth) else auth
                response = self.app.delete(test_url, headers=headers)
                self.assertEqual(expected_status_code, response.status_code)


class TestGetDatasets(BaseAuthAPITest):
    def test_get_dataset_in_a_collection_200(self):
        dataset = self.generate_dataset(name="test")
        test_url = f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}"

        response = self.app.get(test_url)
        print(response)
        self.assertEqual(200, response.status_code)
        self.assertEqual(dataset.dataset_id, response.json["id"])

    def test_get_dataset_shape(self):
        dataset = self.generate_dataset(name="test")
        test_url = f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}"
        response = self.app.get(test_url)
        print(response.json)
        self.assertEqual("test", response.json["title"])

    def test_get_dataset_is_primary_data_shape(self):
        tests = [
            ("PRIMARY", [True]),
            ("SECONDARY", [False]),
            ("BOTH", [True, False]),
        ]
        for is_primary_data, result in tests:
            with self.subTest(f"{is_primary_data}=={result}"):
                metadata = self.sample_dataset_metadata
                metadata.is_primary_data=is_primary_data
                dataset = self.generate_dataset(metadata=metadata)
                test_url = f"/curation/v1/collections/{dataset.collection_id}/datasets/{dataset.dataset_version_id}"
                response = self.app.get(test_url)
                self.assertEqual(result, response.json["is_primary_data"])

    def test_get_nonexistent_dataset_404(self):
        collection = self.generate_unpublished_collection()
        test_url = f"/curation/v1/collections/{collection.collection_id}/datasets/1234-1243-2134-234-1342"
        response = self.app.get(test_url)
        self.assertEqual(404, response.status_code)

    def test_get_datasets_nonexistent_collection_404(self):
        test_url = "/curation/v1/collections/nonexistent/datasets/1234-1243-2134-234-1342"
        headers = self.make_owner_header()
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(404, response.status_code)


class TestPostDataset(BaseAuthAPITest):

    def test_post_datasets_nonexistent_collection_403(self):
        test_url = "/curation/v1/collections/nonexistent/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)

    def test_post_datasets_201(self):
        collection = self.generate_unpublished_collection()
        test_url = f"/curation/v1/collections/{collection.version_id}/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)
        self.assertTrue(response.json["id"])

    def test_post_datasets_super(self):
        collection = self.generate_unpublished_collection()
        test_url = f"/curation/v1/collections/{collection.version_id}/datasets"
        headers = self.make_super_curator_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)

    def test_post_datasets_not_owner_201(self):
        collection = self.generate_collection()
        test_url = f"/curation/v1/collections/{collection.version_id}/datasets"
        headers = self.make_not_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)

    def test_post_datasets_public_collection_405(self):
        collection = self.generate_collection(visibility="PUBLIC")
        test_url = f"/curation/v1/collections/{collection.version_id}/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(405, response.status_code)

    def test_post_datasets_no_auth_401(self):
        collection = self.generate_collection(visibility="PUBLIC")
        test_url = f"/curation/v1/collections/{collection.version_id}/datasets"
        response = self.app.post(test_url)
        self.assertEqual(401, response.status_code)


if __name__ == "__main__":
    unittest.main()
