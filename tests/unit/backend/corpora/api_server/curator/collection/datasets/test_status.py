from backend.corpora.common.corpora_orm import CollectionVisibility, ProcessingStatus
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestDatasetStatus(BaseAuthAPITest):
    def setUp(self):
        super().setUp()
        self.curator_tag = "tag.h5ad"
        self.processing_status = {
            "processing_status": ProcessingStatus.PENDING.name,
        }
        self.queries = ["dataset_uuid", "curator_tag"]
        self.auth_credentials = [
            (self.make_super_curator_header, "super"),
            (self.get_auth_headers, "owner"),
            (None, "none"),
            (self.make_not_owner_header, "not_owner"),
        ]
        self.collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        self.test_url = f"/curation/v1/collections/{self.collection.id}/datasets/status"

    def test__get_status__200(self):
        """Return 200 when the dataset exists, and is part of the collection"""
        dataset = self.generate_dataset(
            self.session,
            collection=self.collection,
            processing_status=self.processing_status,
            curator_tag=self.curator_tag,
        )
        for query in self.queries:
            for auth, auth_description in self.auth_credentials:
                with self.subTest(f"{query} {auth_description} 200"):
                    headers = auth() if callable(auth) else auth
                    query_string = {query: dataset.id} if query == "dataset_uuid" else {query: self.curator_tag}
                    response = self.app.get(self.test_url, headers=headers, query_string=query_string)
                    self.assertEqual(200, response.status_code)
                    self.assertEqual(self.processing_status, response.json)

    def test__get_status__404(self):
        """Return 404 when the dataset does not exist."""
        for query in self.queries:
            for auth, auth_description in self.auth_credentials:
                with self.subTest(f"{query} {auth_description} 404"):
                    headers = auth() if callable(auth) else auth
                    query_string = {query: "fake_id"} if query == "dataset_uuid" else {query: "fake_tag"}
                    response = self.app.get(self.test_url, headers=headers, query_string=query_string)
                    self.assertEqual(404, response.status_code)

    def test__get_status_wrong_collection__404(self):
        """Return 404 when the dataset exists but is not part of the collection."""
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset(
            self.session,
            collection=collection,
            processing_status=self.processing_status,
            curator_tag=self.curator_tag,
        )
        for query in self.queries:
            for auth, auth_description in self.auth_credentials:
                with self.subTest(f"{query} {auth_description} 404"):
                    headers = auth() if callable(auth) else auth
                    query_string = {query: dataset.id} if query == "dataset_uuid" else {query: self.curator_tag}
                    response = self.app.get(self.test_url, headers=headers, query_string=query_string)
                    self.assertEqual(404, response.status_code)

    def test__get_status__400(self):
        """Return 400 when the query parameters are not provided."""
        for query in self.queries:
            for auth, auth_description in self.auth_credentials:
                with self.subTest(f"{query} {auth_description} 400"):
                    headers = auth() if callable(auth) else auth
                    response = self.app.get(self.test_url, headers=headers)
                    self.assertEqual(400, response.status_code)
