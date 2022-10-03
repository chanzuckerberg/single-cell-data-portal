import unittest
from backend.corpora.common.corpora_orm import CollectionVisibility, UploadStatus, IsPrimaryData
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestDeleteDataset(BaseAuthAPITest):
    def test__delete_dataset(self):
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        auth_credentials = [
            (self.make_super_curator_header, "super", 202),
            (self.make_owner_header, "owner", 202),
            (None, "none", 401),
            (self.make_not_owner_header, "not_owner", 403),
        ]
        for auth, auth_description, expected_status_code in auth_credentials:
            with self.subTest(f"{auth_description} {expected_status_code}"):
                collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
                dataset = self.generate_dataset(
                    self.session,
                    collection=collection,
                    processing_status=processing_status,
                )
                test_url = f"/curation/v1/collections/{collection.id}/datasets/{dataset.id}"
                headers = auth() if callable(auth) else auth
                response = self.app.delete(test_url, headers=headers)
                self.assertEqual(expected_status_code, response.status_code)


class TestGetDatasets(BaseAuthAPITest):
    def test_get_dataset_in_a_collection_200(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(self.session, collection=collection, name="test")
        test_url = f"/curation/v1/collections/{collection.id}/datasets/{dataset.id}"

        response = self.app.get(test_url)
        self.assertEqual(200, response.status_code)
        self.assertEqual(dataset.id, response.json["id"])

    def test_get_dataset_shape(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(self.session, collection=collection, name="test")
        test_url = f"/curation/v1/collections/{collection.id}/datasets/{dataset.id}"
        response = self.app.get(test_url)
        self.assertEqual(dataset.name, response.json["title"])

    def test_get_dataset_is_primary_data_shape(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        tests = [
            (IsPrimaryData.PRIMARY, [True]),
            (IsPrimaryData.SECONDARY, [False]),
            (IsPrimaryData.BOTH, [True, False]),
        ]
        for is_primary_data, result in tests:
            with self.subTest(f"{is_primary_data}=={result}"):
                dataset = self.generate_dataset(self.session, collection=collection, is_primary_data=is_primary_data)
                test_url = f"/curation/v1/collections/{collection.id}/datasets/{dataset.id}"
                response = self.app.get(test_url)
                self.assertEqual(result, response.json["is_primary_data"])

    def test_get_nonexistent_dataset_404(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        test_url = f"/curation/v1/collections/{collection.id}/datasets/1234-1243-2134-234-1342"
        response = self.app.get(test_url)
        self.assertEqual(404, response.status_code)

    def test_get_tombstoned_dataset_in_a_collection_404(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(self.session, collection=collection, tombstone=True)
        test_url = f"/curation/v1/collections/{collection.id}/datasets/{dataset.id}"
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
        collection = self.generate_collection(self.session)
        test_url = f"/curation/v1/collections/{collection.id}/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)
        self.assertTrue(response.json["id"])

    def test_post_datasets_super(self):
        collection = self.generate_collection(self.session)
        test_url = f"/curation/v1/collections/{collection.id}/datasets"
        headers = self.make_super_curator_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)

    def test_post_datasets_not_owner_201(self):
        collection = self.generate_collection(self.session)
        test_url = f"/curation/v1/collections/{collection.id}/datasets"
        headers = self.make_not_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)

    def test_post_datasets_public_collection_405(self):
        collection = self.generate_collection(self.session, visibility="PUBLIC")
        test_url = f"/curation/v1/collections/{collection.id}/datasets"
        headers = self.make_owner_header()
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(405, response.status_code)

    def test_post_datasets_no_auth_401(self):
        collection = self.generate_collection(self.session, visibility="PUBLIC")
        test_url = f"/curation/v1/collections/{collection.id}/datasets"
        response = self.app.post(test_url)
        self.assertEqual(401, response.status_code)


if __name__ == "__main__":
    unittest.main()
