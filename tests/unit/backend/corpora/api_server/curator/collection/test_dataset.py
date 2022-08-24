import unittest
from backend.corpora.common.corpora_orm import CollectionVisibility, UploadStatus, IsPrimaryData
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestDeleteDataset(BaseAuthAPITest):
    def test__delete_dataset(self):
        curator_tag = "tag.h5ad"
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        queries = ["dataset_id", "curator_tag"]
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
                        curator_tag=curator_tag,
                    )
                    test_url = f"/curation/v1/collections/{collection.id}/datasets"
                    if query == "dataset_id":
                        query_string = {query: dataset.id}
                    elif query == "curator_tag":
                        query_string = {query: curator_tag}
                    headers = auth() if callable(auth) else auth
                    response = self.app.delete(test_url, headers=headers, query_string=query_string)
                    self.assertEqual(expected_status_code, response.status_code)


class TestPatchDataset(BaseAuthAPITest):
    def test__patch_dataset(self):
        new_tag = "new_tag.h5ad"
        starting_curator_tags = ["tag.h5ad", None]
        test_body = {"curator_tag": new_tag}
        queries = ["dataset_id", "curator_tag"]
        auth_credentials = [
            (self.make_super_curator_header, "super", 204),
            (self.make_owner_header, "owner", 204),
            (None, "none", 401),
            (self.make_not_owner_header, "not_owner", 403),
        ]
        for query in queries:
            for starting_tag in starting_curator_tags:
                if not starting_tag and query == "curator_tag":
                    continue
                for auth, auth_description, expected_status_code in auth_credentials:
                    with self.subTest(f"{query=}, {starting_tag=}, auth={auth_description}, {expected_status_code=}"):
                        collection = self.generate_collection(
                            self.session, visibility=CollectionVisibility.PRIVATE.name
                        )
                        dataset = self.generate_dataset(
                            self.session,
                            collection=collection,
                            curator_tag=starting_tag,
                        )
                        test_url = f"/curation/v1/collections/{collection.id}/datasets"
                        if query == "dataset_id":
                            query_string = {query: dataset.id}
                        elif query == "curator_tag":
                            query_string = {query: starting_tag}
                        headers = auth() if callable(auth) else auth
                        response = self.app.patch(test_url, headers=headers, query_string=query_string, json=test_body)
                        self.assertEqual(expected_status_code, response.status_code)
                        self.session.expire_all()
                        if expected_status_code == 204:
                            self.assertEqual(dataset.curator_tag, new_tag)
                        else:
                            self.assertEqual(dataset.curator_tag, starting_tag)

    def test___invalid_curator_tag(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        test_url = f"/curation/v1/collections/{collection.id}/datasets"

        def _test(_tag, _dataset):
            test_body = {"curator_tag": _tag}
            query_string = {"dataset_id": _dataset.id}
            headers = self.make_owner_header()
            response = self.app.patch(test_url, headers=headers, query_string=query_string, json=test_body)
            self.assertEqual(400, response.status_code)
            self.session.expire_all()
            self.assertIsNone(_dataset.curator_tag)

        tests = [("new", "bad"), ("id", "h5ad"), ("id", "bad")]
        for tag_name, tag_extension in tests:
            with self.subTest(tag_name):
                dataset = self.generate_dataset(self.session, collection=collection)
                _test(tag_name, dataset)
            with self.subTest([tag_name, tag_extension]):
                dataset = self.generate_dataset(self.session, collection=collection)
                name = dataset.id if tag_name == "id" else tag_extension
                tag_name = f"{name}.{tag_extension}"
                _test(tag_name, dataset)

    def test___conflict_curator_tag(self):
        test_tag = "tag.h5ad"
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset_1 = self.generate_dataset(self.session, collection=collection, curator_tag=test_tag)
        dataset_2 = self.generate_dataset(self.session, collection=collection)
        test_body = {"curator_tag": test_tag}
        test_url = f"/curation/v1/collections/{collection.id}/datasets"
        query_string = {"dataset_id": dataset_2.id}
        headers = self.make_owner_header()
        response = self.app.patch(test_url, headers=headers, query_string=query_string, json=test_body)
        self.assertEqual(409, response.status_code)
        self.session.expire_all()
        self.assertEqual(test_tag, dataset_1.curator_tag)
        self.assertIsNone(dataset_2.curator_tag)


class TestGetDatasets(BaseAuthAPITest):
    def test_get_dataset_in_a_collection_200(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(self.session, collection=collection, curator_tag="tag.h5ad", name="test")
        test_url = f"/curation/v1/collections/{collection.id}/datasets"

        test_query_strings = [{"dataset_id": dataset.id}, {"curator_tag": dataset.curator_tag}]
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

        query_string = {"curator_tag": "nonexistent.h5ad"}
        response = self.app.get(test_url, query_string=query_string)
        self.assertEqual(404, response.status_code)

    def test_get_tombstoned_dataset_in_a_collection_404(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(self.session, collection=collection, curator_tag="tag.h5ad", tombstone=True)
        test_url = f"/curation/v1/collections/{collection.id}/datasets"

        test_query_strings = [{"dataset_id": dataset.id}, {"curator_tag": dataset.curator_tag}]
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
