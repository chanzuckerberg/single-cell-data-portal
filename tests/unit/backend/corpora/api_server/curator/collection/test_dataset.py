import unittest

from backend.corpora.common.corpora_orm import CollectionVisibility, UploadStatus
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestDeleteDataset(BaseAuthAPITest):
    def test__delete_dataset(self):
        curator_tag = "tag.h5ad"
        processing_status = {"upload_status": UploadStatus.UPLOADING, "upload_progress": 10.0}
        queries = ["dataset_uuid", "curator_tag"]
        auth_credentials = [
            (self.make_super_curator_header, "super", 202),
            (self.get_auth_headers, "owner", 202),
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
                    test_url = f"/curation/v1/collections/{collection.id}/datasets/"
                    if query == "dataset_uuid":
                        query_string = {query: dataset.id}
                    elif query == "curator_tag":
                        query_string = {query: curator_tag}
                    headers = auth() if callable(auth) else auth
                    response = self.app.delete(test_url, headers=headers, query_string=query_string)
                    self.assertEqual(expected_status_code, response.status_code)


if __name__ == "__main__":
    unittest.main()
