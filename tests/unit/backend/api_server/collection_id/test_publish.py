import json
from datetime import datetime
from unittest.mock import Mock, patch

from backend.common.corpora_orm import CollectionLinkType
from backend.common.entities import Collection
from unit.backend.api_server.base_api_test import BaseAuthAPITest, get_cxguser_token


class TestPublish(BaseAuthAPITest):
    """Test case for publishing a collection."""

    def setUp(self):
        super().setUp()
        self.base_path = "/dp/v1/collections"
        self.headers_authed = {
            "host": "localhost",
            "Content-Type": "application/json",
            "Cookie": get_cxguser_token(),
        }
        self.headers_unauthed = {"host": "localhost", "Content-Type": "application/json"}
        self.mock_published_at = datetime(2000, 12, 25, 0, 0)

    def verify_publish_collection(
        self, pub_collection_id: str, id_to_publish: str = None, mock_timestamp: datetime = None
    ) -> dict:
        """
        Verify publish collection.
        :return: Jsonified response of GET collection/<collection_id>.
        """
        if not mock_timestamp:
            mock_timestamp = self.mock_published_at

        id_to_publish = id_to_publish if id_to_publish else pub_collection_id

        # Publish collection
        body = {"data_submission_policy_version": "1.0"}
        path = f"{self.base_path}/{id_to_publish}/publish"
        with patch("backend.common.entities.collection.datetime") as mock_dt:
            mock_dt.utcnow = Mock(return_value=mock_timestamp)
            response = self.app.post(path, headers=self.headers_authed, data=json.dumps(body))
        self.assertEqual(202, response.status_code)

        self.assertDictEqual({"collection_id": pub_collection_id, "visibility": "PUBLIC"}, json.loads(response.data))
        self.addCleanup(self.delete_collection, pub_collection_id)

        # Cannot call publish for an already published collection
        response = self.app.post(path, headers=self.headers_authed, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

        # Check that the published collection is listed in /collections
        response = self.app.get(self.base_path, headers=self.headers_unauthed)
        self.assertEqual(200, response.status_code)

        ids = [col["id"] for col in json.loads(response.data)["collections"]]
        self.assertIn(pub_collection_id, ids)

        # Check GET collection/<collection_id>
        path = f"{self.base_path}/{pub_collection_id}"
        response = self.app.get(path, headers=self.headers_unauthed)
        self.assertEqual(200, response.status_code)

        response_json = json.loads(response.data)
        self.assertEqual("PUBLIC", response_json["visibility"])
        self.assertEqual(pub_collection_id, response_json["id"])

        return response_json

    def test__publish_a_new_collection__OK(self):
        """Publish a new collection with a single dataset."""
        collection = self.generate_collection(self.session)
        self.generate_dataset(self.session, collection_id=collection.id)

        self.assertIsNone(collection.published_at)
        response_json = self.verify_publish_collection(collection.id, mock_timestamp=self.mock_published_at)

        # Check collection published_at
        self.assertEqual(self.mock_published_at, datetime.utcfromtimestamp(response_json["published_at"]))

        dataset = response_json["datasets"][0]
        self.assertEqual(self.mock_published_at, datetime.utcfromtimestamp(dataset["published_at"]))

    def test__published_at_not_updated_for_publish_a_collection_revision(self):
        """Publish a collection revision should not update published_at timestamp."""
        # Publish a new collection
        collection = self.generate_collection(self.session)
        self.generate_dataset(self.session, collection_id=collection.id)

        self.assertIsNone(collection.published_at)
        response_json = self.verify_publish_collection(collection.id, mock_timestamp=self.mock_published_at)
        self.assertEqual(self.mock_published_at, datetime.utcfromtimestamp(response_json["published_at"]))

        # Start a revision of the collection
        response = self.app.post(f"{self.base_path}/{collection.id}", headers=self.headers_authed)
        self.assertEqual(201, response.status_code)

        response_json = json.loads(response.data)
        revision_id = response_json["id"]
        self.assertEqual("PRIVATE", response_json["visibility"])

        # Publish revision
        mock_revision_published_dt = datetime(2001, 1, 25, 0, 0)
        response_json = self.verify_publish_collection(
            collection.id, id_to_publish=revision_id, mock_timestamp=mock_revision_published_dt
        )

        # published_at should not be updated - only updates on initial publish
        self.assertEqual(self.mock_published_at, datetime.utcfromtimestamp(response_json["published_at"]))

        dataset = response_json["datasets"][0]
        self.assertEqual(self.mock_published_at, datetime.utcfromtimestamp(dataset["published_at"]))

    def test__publish_collection_w_multiple_datasets__OK(self):
        """Publish collection with multiple datasets."""
        collection_id = self.generate_collection(self.session).id
        dataset_ids = [
            self.generate_dataset(self.session, collection_id=collection_id).id,
            self.generate_dataset(self.session, collection_id=collection_id).id,
        ]

        response_json = self.verify_publish_collection(collection_id)
        self.assertEqual(self.mock_published_at, datetime.utcfromtimestamp(response_json["published_at"]))

        # Check datasets
        res_datasets = response_json["datasets"]
        res_dataset_ids = [dataset["id"] for dataset in res_datasets]

        self.assertListEqual(sorted(dataset_ids), sorted(res_dataset_ids))
        for dataset in res_datasets:
            self.assertTrue(dataset["published"])
            self.assertEqual(self.mock_published_at, datetime.utcfromtimestamp(dataset["published_at"]))

    def verify_publish_collection_with_links(self, collection: Collection, id_to_publish: str = None):
        link_names = [link.link_name for link in collection.links]

        response_json = self.verify_publish_collection(collection.id, id_to_publish=id_to_publish)
        self.assertEqual(self.mock_published_at, datetime.utcfromtimestamp(response_json["published_at"]))

        # Check links
        res_links = [link["link_name"] for link in response_json["links"]]
        self.assertListEqual(sorted(link_names), sorted(res_links))

    def test__publish_collection_with_links__OK(self):
        """Publish collection with a link."""
        collection = self.generate_collection(
            self.session,
            links=[
                {"link_name": "test_link", "link_type": CollectionLinkType.PROTOCOL, "link_url": "https://link.link"}
            ],
            published_at=self.mock_published_at,
        )
        self.generate_dataset(self.session, collection_id=collection.id, published_at=self.mock_published_at)

        self.verify_publish_collection_with_links(collection)

    def test__publish_collection_revision_with_links__OK(self):
        collection = Collection.get_collection(self.session, collection_id="test_collection_with_link")
        self.generate_dataset(self.session, collection_id=collection.id, published_at=self.mock_published_at)
        revision = collection.create_revision()

        collection.update(published_at=self.mock_published_at, keep_links=True)

        self.verify_publish_collection_with_links(collection, revision.id)

    def test__publish_collection_revision_with_links_and_dataset_changes__OK(self):
        collection = Collection.get_collection(
            self.session, collection_id="test_collection_with_link_and_dataset_changes"
        )
        self.generate_dataset(self.session, collection_id=collection.id, published_at=self.mock_published_at)
        revision = collection.create_revision()
        self.generate_dataset(self.session, collection_id=revision.id)  # Collection will have Dataset changes

        collection.update(published_at=self.mock_published_at, keep_links=True)

        self.verify_publish_collection_with_links(collection, revision.id)

    @patch("backend.common.utils.cloudfront.create_invalidation_for_index_paths")
    def test_publish_collection_does_cloudfront_invalidation(self, mock_cloudfront):
        """Publish a new collection with a single dataset."""
        collection = self.generate_collection(self.session)
        self.generate_dataset(self.session, collection_id=collection.id)

        self.assertIsNone(collection.published_at)
        self.verify_publish_collection(collection.id, mock_timestamp=self.mock_published_at)

        mock_cloudfront.assert_called_once()

    def test__not_owner__403(self):
        """Publish a collection as a non-owner."""
        collection_id = self.generate_collection(self.session, owner="someone_else").id
        path = f"{self.base_path}/{collection_id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        response = self.app.post(path, headers=self.headers_authed, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    def test__bad_id__403(self):
        """Publish a collection with a bad uuid."""
        collection_id = "bad_id"
        body = {"data_submission_policy_version": "1.0"}
        path = f"{self.base_path}/{collection_id}/publish"
        response = self.app.post(path, headers=self.headers_authed, data=json.dumps(body))
        self.assertEqual(403, response.status_code)


class TestPublishCurators(BaseAuthAPITest):
    def verify_can_publish(self, collection, user):
        self.generate_dataset(self.session, collection_id=collection.id)
        path = f"/dp/v1/collections/{collection.id}/publish"
        body = {"data_submission_policy_version": "1.0"}
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_cxguser_token(user)}
        response = self.app.post(path, headers=headers, data=json.dumps(body))
        self.assertEqual(202, response.status_code)

    def test__can_publish_owned_collection(self):
        collection = self.generate_collection(self.session)
        self.verify_can_publish(collection, "owner")

    def test__can_publish_non_owned_collection_as_super_curator(self):
        collection = self.generate_collection(self.session, owner="someone_else")
        self.verify_can_publish(collection, "super")
