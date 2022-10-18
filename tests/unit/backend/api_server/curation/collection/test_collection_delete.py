import unittest

from backend.common.corpora_orm import CollectionVisibility
from unit.backend.api_server.base_api_test import BaseAuthAPITest


class TestDeleteCollection(BaseAuthAPITest):
    def _test(self, collection_id, header, expected_status):
        if header == "owner":
            headers = self.make_owner_header()
        elif header == "super":
            headers = self.make_super_curator_header()
        elif header == "not_owner":
            headers = self.make_not_owner_header()
        elif "noauth":
            headers = {}

        response = self.app.delete(f"/curation/v1/collections/{collection_id}", headers=headers)
        self.assertEqual(expected_status, response.status_code)
        if response.status_code == 204:
            response = self.app.delete(f"/curation/v1/collections/{collection_id}", headers=headers)
            self.assertEqual(403, response.status_code)

    def test__delete_public_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 405), ("super", 405)]
        public_collection_id = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name).id
        for auth, expected_response in tests:
            with self.subTest(auth):
                self._test(public_collection_id, auth, expected_response)

    def test__delete_revision_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 204), ("super", 204)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                public_collection_id = self.generate_collection(
                    self.session, visibility=CollectionVisibility.PUBLIC.name
                ).id
                revision_collection_id = self.generate_collection(
                    self.session, visibility=CollectionVisibility.PRIVATE.name, revision_of=public_collection_id
                ).id
                self._test(revision_collection_id, auth, expected_response)

    def test__delete_private_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 204), ("super", 204)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                private_collection_id = self.generate_collection(
                    self.session, visibility=CollectionVisibility.PRIVATE.name
                ).id
                self._test(private_collection_id, auth, expected_response)

    def test__delete_tombstone_collection(self):
        tests = [("not_owner", 403), ("noauth", 401), ("owner", 403), ("super", 403)]
        for auth, expected_response in tests:
            with self.subTest(auth):
                tombstone_collection_id = self.generate_collection(
                    self.session, visibility=CollectionVisibility.PUBLIC.name, tombstone=True
                ).id
                self._test(tombstone_collection_id, auth, expected_response)


if __name__ == "__main__":
    unittest.main()
