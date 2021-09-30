import json

from backend.corpora.common.corpora_orm import CollectionVisibility, CollectionLinkType
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest, BasicAuthAPITestCurator
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token


class TestPublish(BaseAuthAPITest):
    """Test case for publishing a collection."""

    def setUp(self):
        super().setUp()
        self.base_path = "/dp/v1/collections"
        self.headers_authed = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": get_auth_token(self.app),
        }
        self.headers_unauthed = {"host": "localhost", "Content-Type": "application/json"}

    def verify_publish_collection(self, collection_id: str) -> dict:
        """
        Verify publish collection.
        Returns:
            response_json (dict): Jsonified response of GET collection/<collection_id>.
        """

        # Publish collection
        path = f"{self.base_path}/{collection_id}/publish"
        response = self.app.post(path, headers=self.headers_authed)
        self.assertEqual(202, response.status_code)

        self.assertDictEqual({"collection_uuid": collection_id, "visibility": "PUBLIC"}, json.loads(response.data))
        self.addCleanup(self.delete_collection, collection_id, "PUBLIC")

        # Cannot call publish for an already published collection
        response = self.app.post(path, headers=self.headers_authed)
        self.assertEqual(403, response.status_code)

        # Check that the published collection is listed in /collections
        response = self.app.get(self.base_path, headers=self.headers_unauthed)
        self.assertEqual(200, response.status_code)

        ids = [col["id"] for col in json.loads(response.data)["collections"]]
        self.assertIn(collection_id, ids)

        # Check GET collection/<collection_id>
        path = f"{self.base_path}/{collection_id}"
        response = self.app.get(path, headers=self.headers_unauthed)
        self.assertEqual(200, response.status_code)

        response_json = json.loads(response.data)
        self.assertEqual("PUBLIC", response_json["visibility"])
        self.assertEqual(collection_id, response_json["id"])

        return response_json

    def test__publish_collection__OK(self):
        """Publish collection with a single dataset."""
        collection = self.generate_collection(
            self.session,
        )
        self.generate_dataset(self.session, collection_id=collection.id, collection_visibility=collection.visibility)
        self.verify_publish_collection(collection.id)

    def test__publish_collection_w_multiple_datasets__OK(self):
        """Publish collection with multiple datasets."""
        collection_id = self.generate_collection(self.session).id
        dataset_ids = [
            self.generate_dataset(
                self.session, collection_id=collection_id, collection_visibility=CollectionVisibility.PRIVATE.name
            ).id,
            self.generate_dataset(
                self.session, collection_id=collection_id, collection_visibility=CollectionVisibility.PRIVATE.name
            ).id,
        ]

        response_json = self.verify_publish_collection(collection_id)

        # Check datasets
        res_datasets = response_json["datasets"]
        res_dataset_ids = [dataset["id"] for dataset in res_datasets]

        self.assertListEqual(sorted(dataset_ids), sorted(res_dataset_ids))
        for dataset in res_datasets:
            self.assertTrue(dataset["published"])

    def test__publish_collection_with_links__OK(self):
        """Publish collection with a link."""
        collection = self.generate_collection(
            self.session,
            links=[
                {"link_name": "test_link", "link_type": CollectionLinkType.PROTOCOL, "link_url": "https://link.link"}
            ],
        )
        link_names = [link.link_name for link in collection.links]
        dataset_id = self.generate_dataset(
            self.session, collection_id=collection.id, collection_visibility=CollectionVisibility.PRIVATE.name
        ).id

        response_json = self.verify_publish_collection(collection.id)
        self.assertEqual(dataset_id, response_json["datasets"][0]["id"])

        # Check links
        res_links = [link["link_name"] for link in response_json["links"]]
        self.assertListEqual(sorted(link_names), sorted(res_links))

    def test__not_owner__403(self):
        """Publish a collection as a non-owner."""
        collection_id = self.generate_collection(self.session, owner="someone_else").id
        path = f"{self.base_path}/{collection_id}/publish"
        response = self.app.post(path, headers=self.headers_authed)
        self.assertEqual(403, response.status_code)

    def test__bad_uuid__403(self):
        """Publish a collection with a bad uuid."""
        collection_id = "bad_uuid"
        path = f"{self.base_path}/{collection_id}/publish"
        response = self.app.post(path, headers=self.headers_authed)
        self.assertEqual(403, response.status_code)


class TestPublishCurators(BasicAuthAPITestCurator):
    def test__can_publish_owned_collection(self):
        collection_id = self.generate_collection(self.session).id
        self.generate_dataset(self.session, collection_id=collection_id, collection_visibility="PRIVATE")
        path = f"/dp/v1/collections/{collection_id}/publish"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.post(path, headers=headers)
        self.assertEqual(202, response.status_code)

    def test__can_publish_non_owned_collection(self):
        collection_id = self.generate_collection(self.session, owner="someone_else").id
        self.generate_dataset(self.session, collection_id=collection_id, collection_visibility="PRIVATE")
        path = f"/dp/v1/collections/{collection_id}/publish"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.post(path, headers=headers)
        self.assertEqual(202, response.status_code)
