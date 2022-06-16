from backend.corpora.common.corpora_orm import CollectionVisibility, ProcessingStatus
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestDatasetStatus(BaseAuthAPITest):
    def test_get_status(self):
        curator_tag = "tag.h5ad"
        processing_status = {
            "processing_status": ProcessingStatus.PENDING.name,
        }
        queries = ["dataset_uuid", "curator_tag"]
        auth_credentials = [
            (self.make_super_curator_header, "super", 200),
            (self.get_auth_headers, "owner", 200),
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
                    test_url = f"/curation/v1/collections/{collection.id}/datasets/status"
                    query_string = {query: dataset.id} if query == "dataset_uuid" else {query: curator_tag}
                    headers = auth() if callable(auth) else auth
                    response = self.app.get(test_url, headers=headers, query_string=query_string)
                    self.assertEqual(expected_status_code, response.status_code)
                    if response.status_code == 200:
                        self.assertEqual(processing_status, response.json)
