import unittest

from backend.corpora.common.corpora_orm import CollectionVisibility
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestDeleteCollection(BaseAuthAPITest):
    # cannot delete public
    # cannot delete a revision
    # can delete a private collection
    def setUp(self):
        super().setUp()
        self.public_collection_uuid = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name
        ).id
        self.revision_collection_uuid = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, revision_of=self.public_collection_uuid
        ).id

    def _test(self, collection_uuid, header, expected_status):
        if header == "owner":
            headers = self.get_auth_headers()
        elif header == "super":
            headers = self.make_super_curator_header()
        elif header == "not_owner":
            headers = self.make_not_owner_header()
        elif "noauth":
            headers = {}

        response = self.app.delete(f"/curation/v1/collections/{collection_uuid}", headers=headers)
        self.assertEqual(expected_status, response.status_code)
        if response.status_code == 204:
            # After deleting the collection is not longer found. No tombstone is created because private collections
            # are deleted, not tombstoned.
            response = self.app.delete(f"/curation/v1/collections/{collection_uuid}", headers=headers)
            self.assertEqual(403, response.status_code)

    def test__delete_public_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 403), ("super", 403)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                self._test(self.public_collection_uuid, auth, expected_response)

    def test__delete_revision_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 405), ("super", 405)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                self._test(self.revision_collection_uuid, auth, expected_response)

    def test__delete_private_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 204), ("super", 204)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                private_collection_uuid = self.generate_collection(
                    self.session, visibility=CollectionVisibility.PRIVATE.name
                ).id
                self._test(private_collection_uuid, auth, expected_response)


if __name__ == "__main__":
    unittest.main()
