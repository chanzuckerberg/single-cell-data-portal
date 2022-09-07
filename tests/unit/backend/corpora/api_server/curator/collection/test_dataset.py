import unittest
from backend.corpora.common.corpora_orm import CollectionVisibility, UploadStatus, IsPrimaryData
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestDeleteDataset(BaseAuthAPITest):
    def test__delete_dataset(self):
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        queries = ["dataset_id"]
        auth_credentials = [
            (self.make_super_curator_header, "super", 202),
            (self.make_owner_header, "owner", 202),
            (None, "none", 401),
            (self.make_not_owner_header, "not_owner", 403),
        ]
        for query in queries:
            for auth, auth_description, expected_status_code in auth_credentials:
                with self.subTest(f"{query} {auth_description} {expected_status_code}"):
                    collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
                    dataset = self.generate_dataset(
                        self.session,
                        collection=collection,
                        processing_status=processing_status,
                    )
                    test_url = f"/curation/v1/collections/{collection.id}/datasets"
                    if query == "dataset_id":
                        query_string = {query: dataset.id}
                    headers = auth() if callable(auth) else auth
                    response = self.app.delete(test_url, headers=headers, query_string=query_string)
                    self.assertEqual(expected_status_code, response.status_code)


class TestGetDatasets(BaseAuthAPITest):
    def test_get_dataset_in_a_collection_200(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(self.session, collection=collection, name="test")
        test_url = f"/curation/v1/collections/{collection.id}/datasets"

        test_query_strings = [{"dataset_id": dataset.id}]
        for query_string in test_query_strings:
            with self.subTest(query_string):
                response = self.app.get(test_url, query_string=query_string)
                self.assertEqual(200, response.status_code)
                self.assertEqual(dataset.id, response.json["id"])

    def test_get_dataset_shape(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        test_url = f"/curation/v1/collections/{collection.id}/datasets"

        with self.subTest("name to title"):
            dataset = self.generate_dataset(self.session, collection=collection, name="test")
            response = self.app.get(test_url, query_string={"dataset_id": dataset.id})
            self.assertEqual(dataset.name, response.json["title"])

    def test_get_dataset_is_primary_data_shape(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        test_url = f"/curation/v1/collections/{collection.id}/datasets"

        tests = [
            (IsPrimaryData.PRIMARY, [True]),
            (IsPrimaryData.SECONDARY, [False]),
            (IsPrimaryData.BOTH, [True, False]),
        ]
        for is_primary_data, result in tests:
            with self.subTest(f"{is_primary_data}=={result}"):
                dataset = self.generate_dataset(self.session, collection=collection, is_primary_data=is_primary_data)
                response = self.app.get(test_url, query_string={"dataset_id": dataset.id})
                self.assertEqual(result, response.json["is_primary_data"])

    def test_get_nonexistent_dataset_404(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        test_url = f"/curation/v1/collections/{collection.id}/datasets"

        query_string = {"dataset_id": "1234-1243-2134-234-1342"}
        response = self.app.get(test_url, query_string=query_string)
        self.assertEqual(404, response.status_code)

    def test_get_tombstoned_dataset_in_a_collection_404(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(self.session, collection=collection, tombstone=True)
        test_url = f"/curation/v1/collections/{collection.id}/datasets"

        test_query_strings = [{"dataset_id": dataset.id}]
        for query_string in test_query_strings:
            with self.subTest(query_string):
                response = self.app.get(test_url, query_string=query_string)
                self.assertEqual(404, response.status_code)

    def test_get_datasets_nonexistent_collection_404(self):
        test_url = "/curation/v1/collections/nonexistent/datasets"
        headers = self.make_owner_header()
        response = self.app.get(test_url, headers=headers)
        self.assertEqual(404, response.status_code)


if __name__ == "__main__":
    unittest.main()
