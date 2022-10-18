from backend.common.corpora_orm import CollectionVisibility
from unit.backend.api_server.base_api_test import BaseAuthAPITest


class TestPostRevision(BaseAuthAPITest):
    def test__post_revision__no_auth(self):
        collection_id = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name).id
        response = self.app.post(f"/curation/v1/collections/{collection_id}/revision")
        self.assertEqual(401, response.status_code)

    def test__post_revision__Not_Public(self):
        collection_id = self.generate_collection(self.session).id
        headers = self.make_super_curator_header()
        response = self.app.post(f"/curation/v1/collections/{collection_id}/revision", headers=headers)
        self.assertEqual(403, response.status_code)

    def test__post_revision__Not_Owner(self):
        collection_id = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="someone else"
        ).id
        response = self.app.post(
            f"/curation/v1/collections/{collection_id}/revision",
            headers=self.make_owner_header(),
        )
        self.assertEqual(403, response.status_code)

    def test__post_revision__OK(self):
        collection_id = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name).id
        response = self.app.post(
            f"/curation/v1/collections/{collection_id}/revision",
            headers=self.make_owner_header(),
        )
        self.assertEqual(201, response.status_code)
        self.assertNotEqual(collection_id, response.json["id"])

    def test__post_revision__Super_Curator(self):
        collection_id = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name).id
        headers = self.make_super_curator_header()
        response = self.app.post(f"/curation/v1/collections/{collection_id}/revision", headers=headers)
        self.assertEqual(201, response.status_code)
        self.assertNotEqual(collection_id, response.json["id"])
