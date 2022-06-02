from backend.corpora.common.corpora_orm import CollectionVisibility
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest


class TestPostRevision(BaseAuthAPITest):
    def test__post_revision__no_auth(self):
        collection_uuid = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name).id
        response = self.app.post(f"/curation/v1/collections/{collection_uuid}/revision")
        self.assertEqual(401, response.status_code)

    def test__post_revision__Not_Public(self):
        collection_uuid = self.generate_collection(self.session).id
        headers = self.make_super_curator_header()
        response = self.app.post(f"/curation/v1/collections/{collection_uuid}/revision", headers=headers)
        self.assertEqual(403, response.status_code)

    def test__post_revision__Not_Owner(self):
        collection_uuid = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC.name, owner="someone else"
        ).id
        response = self.app.post(
            f"/curation/v1/collections/{collection_uuid}/revision",
            headers=self.get_auth_headers(),
        )
        self.assertEqual(403, response.status_code)

    def test__post_revision__OK(self):
        collection_uuid = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name).id
        response = self.app.post(
            f"/curation/v1/collections/{collection_uuid}/revision",
            headers=self.get_auth_headers(),
        )
        self.assertEqual(201, response.status_code)
        self.assertNotEqual(collection_uuid, response.json["revision_id"])

    def test__post_revision__Super_Curator(self):
        collection_uuid = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name).id
        headers = self.make_super_curator_header()
        response = self.app.post(f"/curation/v1/collections/{collection_uuid}/revision", headers=headers)
        self.assertEqual(201, response.status_code)
        self.assertNotEqual(collection_uuid, response.json["revision_id"])
