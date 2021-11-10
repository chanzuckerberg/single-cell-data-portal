import json
import random
import string
import typing
from datetime import datetime
from mock import Mock, patch
from unittest import mock

from sqlalchemy.exc import SQLAlchemyError

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Dataset, Collection
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.exceptions import CorporaException
from backend.corpora.common.utils.json import CustomJSONEncoder
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class BaseRevisionTest(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        super().setUp()
        pub_collection = self.generate_collection(self.session, visibility="PUBLIC")
        for i in range(2):
            self.generate_dataset_with_s3_resources(
                self.session, collection_visibility="PUBLIC", collection_id=pub_collection.id, published=True
            )
        self.pub_collection = pub_collection
        self.rev_collection = pub_collection.revision()
        self.headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

    def assertPublishedCollectionOK(self, expected_body, s3_objects):
        """Checks that the published collection is as expected and S3 Objects exist"""
        with self.subTest("published artifacts and explorer s3 object ok"):
            for bucket, file_name in s3_objects:
                self.assertS3FileExists(bucket, file_name)
        with self.subTest("publish collection ok"):
            resp = self.app.get(f"/dp/v1/collections/{self.pub_collection.id}")
            self.assertEqual(200, resp.status_code)
            actual_body = json.loads(resp.data)
            for link in expected_body.pop("links"):
                self.assertIn(link, actual_body["links"])
            for keys in expected_body.keys():
                self.assertEqual(expected_body[keys], actual_body[keys])

    def refresh_datasets(self) -> typing.List[Dataset]:
        for dataset in self.rev_collection.datasets:
            Dataset(dataset).delete()
        refreshed_datasets = []
        for dataset in self.pub_collection.datasets:
            ds = self.generate_dataset_with_s3_resources(
                self.session,
                collection_visibility="PRIVATE",
                collection_id=self.rev_collection.id,
                original_id=dataset.id,
                revision=dataset.revision + 1,
            )
            refreshed_datasets.append(ds)
        return refreshed_datasets

    def get_s3_objects_from_collections(self) -> typing.Tuple[typing.List, typing.List]:
        """
        :return: a list of s3 objects in the published collection, and a list of s3 objects the revision collection.
        """
        rev_s3_objects = []
        pub_s3_objects = []

        for i in range(len(self.pub_collection.datasets)):
            pub_s3_objects.extend(self.get_s3_object_paths_from_dataset(self.pub_collection.datasets[i]))
            rev_s3_objects.extend(self.get_s3_object_paths_from_dataset(self.rev_collection.datasets[i]))
        return pub_s3_objects, rev_s3_objects


class TestRevision(BaseRevisionTest):
    """Test case for starting a collection's revision."""

    def verify_start_revision(self, collection_id: str) -> dict:
        """
        Verify start of a collection's revision.
        :return: Jsonified response of POST collection/<collection_id>.
        """
        path = f"/dp/v1/collections/{collection_id}"
        response = self.app.post(path, self.headers)
        self.assertEqual(201, response.status_code)

        response_json = json.loads(response.data)
        self.assertEqual("PRIVATE", response_json["visibility"])

        return response_json

    def verify_get_revision(self, collection_id: str, dataset_ids: typing.List[str] = None) -> dict:
        """
        Verify the contents of a collection under revision.
        :return: Jsonified response of GET collection/<collection_id>?visibility=PRIVATE.
        """
        path = f"/dp/v1/collections/{collection_id}?visibility=PRIVATE"
        response = self.app.get(path, headers=self.headers)
        response_json = json.loads(response.data)

        with self.subTest("Test datasets in revised collection are not original datasets"):
            new_dataset_ids = [x["id"] for x in response_json["datasets"]]
            for dataset_id in dataset_ids:
                self.assertNotIn(dataset_id, new_dataset_ids)

        with self.subTest("Test revised datasets point at original datasets"):
            original_dataset_ids = [x["original_id"] for x in response_json["datasets"]]
            for dataset_id in dataset_ids:
                self.assertIn(dataset_id, original_dataset_ids)

        with self.subTest("Check assets point at revised dataset"):
            for dataset in response_json["datasets"]:
                dataset_id = dataset["id"]
                asset_dataset_ids = {asset["dataset_id"] for asset in dataset["dataset_assets"]}
                self.assertEqual(dataset_id, asset_dataset_ids.pop())

        with self.subTest("Test revised_at not updated under revision"):
            # Collection revised_at
            self.assertIsNone(response_json.get("revised_at"))

            # Dataset revised_at
            for dataset in response_json["datasets"]:
                self.assertIsNone(dataset.get("revised_at"))

        return response_json

    def verify_unauthed_get_revision(self, collection_id: str, expected_body: dict) -> None:
        """Verify unauthorized view of a collection under revision."""

        path = f"/dp/v1/collections/{collection_id}?visibility=PRIVATE"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(path, headers=headers)
        self.assertEqual(200, response.status_code)

        response_json = json.loads(response.data)
        self.assertEqual("READ", response_json.pop("access_type"))
        self.assertEqual(expected_body, response_json)

    def test__start_revision_of_a_collection__201(self):
        """Start a revision of a collection."""
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC)
        dataset_ids = [
            self.generate_dataset_with_s3_resources(self.session, collection=collection, published=True).id,
            self.generate_dataset_with_s3_resources(self.session, collection=collection, published=True).id,
        ]

        res_post_json = self.verify_start_revision(collection.id)
        res_get_json = self.verify_get_revision(collection.id, dataset_ids)
        self.assertEqual(res_post_json, res_get_json)

        self.assertEqual("WRITE", res_get_json.pop("access_type"))
        self.verify_unauthed_get_revision(collection.id, res_get_json)

    def test__start_revision_of_a_collection_w_links__201(self):
        """Start a revision of a collection with links."""
        collection = self.generate_collection(
            self.session,
            visibility=CollectionVisibility.PUBLIC,
            links=[
                {"link_name": "Link 1", "link_url": "This is a new link", "link_type": "OTHER"},
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
            ],
        )
        dataset_ids = [
            self.generate_dataset_with_s3_resources(self.session, collection=collection, published=True).id,
            self.generate_dataset_with_s3_resources(self.session, collection=collection, published=True).id,
        ]
        link_names = [link.link_name for link in collection.links]

        res_post_json = self.verify_start_revision(collection.id)

        # Check link names
        res_links = [link["link_name"] for link in res_post_json["links"]]
        self.assertListEqual(sorted(link_names), sorted(res_links))

        res_get_json = self.verify_get_revision(collection.id, dataset_ids)
        self.assertEqual(res_post_json, res_get_json)

        self.assertEqual("WRITE", res_get_json.pop("access_type"))
        self.verify_unauthed_get_revision(collection.id, res_get_json)

    def test__revision__409(self):
        """Starting a revision on a revision."""
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC)
        test_url = f"/dp/v1/collections/{collection.id}"

        # Start a revision
        response = self.app.post(test_url, headers=self.headers)
        self.assertEqual(201, response.status_code)

        # Try to start a revision again
        response = self.app.post(test_url, headers=self.headers)
        self.assertEqual(409, response.status_code)

    def test__revision_nonexistent__403(self):
        """Start a revision on a non-existing collection."""
        response = self.app.post("/dp/v1/collections/random", headers=self.headers)
        self.assertEqual(403, response.status_code)

    def test__revision_not_owner__403(self):
        """Start a revision on a collection as a non-owner."""
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC, owner="someone else"
        )
        test_url = f"/dp/v1/collections/{collection.id}"

        response = self.app.post(test_url, headers=self.headers)
        self.assertEqual(403, response.status_code)


class TestDeleteRevision(BaseRevisionTest):
    """Test case for deleting a collection or datasets under revision."""

    def setUp(self):
        super().setUp()
        url = f"/dp/v1/collections/{self.pub_collection.id}"
        self.test_url_collect_private = f"{url}?visibility=PRIVATE"
        self.test_url_collection_public = f"{url}?visibility=PUBLIC"

    def test__revision_deleted__204(self):
        """Delete a collection under revision."""
        # Delete the revision
        resp = self.app.delete(self.test_url_collect_private, headers=self.headers)
        self.assertEqual(204, resp.status_code)

        # Cannot get the revision
        resp = self.app.get(self.test_url_collect_private, headers=self.headers)
        self.assertEqual(403, resp.status_code)

    def test__revision_deleted_with_published_datasets(self):
        """
        The published dataset artifacts should be intact after deleting a
        collection revision.
        """
        expected_body = json.loads(json.dumps(self.pub_collection.reshape_for_api(), cls=CustomJSONEncoder))
        pub_s3_objects, rev_s3_objects = self.get_s3_objects_from_collections()

        # Revision and Published collection refer to the same S3 resources
        self.assertEqual(pub_s3_objects, rev_s3_objects)

        # Delete the revision
        resp = self.app.delete(self.test_url_collect_private, headers=self.headers)
        self.assertEqual(204, resp.status_code)
        self.assertPublishedCollectionOK(expected_body, pub_s3_objects)

    def test__revision_deleted_with_new_datasets(self):
        """The new datasets should be deleted when the revison is deleted."""
        # Generate revision dataset
        rev_dataset = self.generate_dataset_with_s3_resources(
            self.session, collection_visibility="PRIVATE", collection_id=self.rev_collection.id, published=False
        )
        s3_objects = self.get_s3_object_paths_from_dataset(rev_dataset)

        # Check resources exist
        with self.subTest("new artifacts and explorer s3 objects exist"):
            for bucket, file_name in s3_objects:
                self.assertS3FileExists(bucket, file_name)

        # Delete Revision
        resp = self.app.delete(self.test_url_collect_private, headers=self.headers)
        self.assertEqual(204, resp.status_code)

        with self.subTest("new artifacts and explorer s3 objects deleted"):
            for bucket, file_name in s3_objects:
                self.assertS3FileDoesNotExist(bucket, file_name)

    def test__revision_deleted_with_refreshed_datasets(self):
        """
        The refreshed datasets should be deleted and the published dataset
        intact. The published dataset artifacts should be intact after
        deleting a collection revision.
        """
        expected_body = json.loads(json.dumps(self.pub_collection.reshape_for_api(), cls=CustomJSONEncoder))
        self.refresh_datasets()
        pub_s3_objects, rev_s3_objects = self.get_s3_objects_from_collections()

        # Refreshed datasets do not point to the published resources in s3.
        for s3_object in rev_s3_objects:
            self.assertNotIn(s3_object, pub_s3_objects)

        # Delete the revision
        resp = self.app.delete(self.test_url_collect_private, headers=self.headers)
        self.assertEqual(204, resp.status_code)

        with self.subTest("refreshed artifacts and explorer s3 objects deleted"):
            for bucket, file_name in rev_s3_objects:
                self.assertS3FileDoesNotExist(bucket, file_name)
        self.assertPublishedCollectionOK(expected_body, pub_s3_objects)

    def test__delete_refreshed_dataset_in_a_revision(self):
        """
        The refreshed datasets should be deleted and the published dataset
        restored in the revision.
        """
        expected_pub_body = json.loads(json.dumps(self.pub_collection.reshape_for_api(), cls=CustomJSONEncoder))
        expected_rev_body = self.remove_timestamps(
            json.loads(json.dumps(self.rev_collection.reshape_for_api(), cls=CustomJSONEncoder)),
            remove=["id", "dataset_id", "dataset_deployments"],
        )
        refreshed_datasets = self.refresh_datasets()
        pub_s3_objects, rev_s3_objects = self.get_s3_objects_from_collections()

        # Refreshed datasets do not point to the published resources in s3.
        for s3_object in rev_s3_objects:
            self.assertNotIn(s3_object, pub_s3_objects)

        # Delete the refreshed_datasets
        for ds in refreshed_datasets:
            resp = self.app.delete(f"/dp/v1/datasets/{ds.id}", headers=self.headers)
            self.assertEqual(202, resp.status_code)

        with self.subTest("refreshed artifacts and explorer s3 objects deleted"):
            for bucket, file_name in rev_s3_objects:
                self.assertS3FileDoesNotExist(bucket, file_name)
        self.assertPublishedCollectionOK(expected_pub_body, pub_s3_objects)

        # Check that the original datasets info has been restored to the revision dataset.
        self.session.expire_all()
        actual_rev_body = self.remove_timestamps(
            json.loads(json.dumps(self.rev_collection.reshape_for_api(), cls=CustomJSONEncoder)),
            remove=["id", "dataset_id", "dataset_deployments"],
        )
        self.assertEqual(expected_rev_body, actual_rev_body)

    def test__delete_published_dataset_during_revision(self):
        """
        The dataset is tombstone in the revision. The published artifacts are
        intact.
        """
        expected_body = json.loads(json.dumps(self.pub_collection.reshape_for_api(), cls=CustomJSONEncoder))
        pub_s3_objects, _ = self.get_s3_objects_from_collections()

        # Delete a published dataset in the revision
        rev_dataset_count = len(self.rev_collection.datasets)
        rev_dataset_id = self.rev_collection.datasets[0].id
        test_dataset_url = f"/dp/v1/datasets/{rev_dataset_id}"
        resp = self.app.delete(test_dataset_url, headers=self.headers)
        self.assertEqual(202, resp.status_code)

        # Get the revision authenticated
        resp = self.app.get(self.test_url_collect_private, headers=self.headers)
        self.assertEqual(200, resp.status_code)
        self.assertEqual(rev_dataset_count, len(json.loads(resp.data)["datasets"]))

        # Get the revision unauthenticated
        resp = self.app.get(self.test_url_collect_private)
        self.assertEqual(200, resp.status_code)
        # The dataset is a tombstone in the revision
        self.assertEqual(rev_dataset_count - 1, len(json.loads(resp.data)["datasets"]))

        self.session.expire_all()
        for dataset in self.rev_collection.datasets:
            if dataset.id == rev_dataset_id:
                self.assertTrue(dataset.tombstone)
                break
        self.assertPublishedCollectionOK(expected_body, pub_s3_objects)


def get_random_etag(*args, **kwargs):
    return {"ETag": ''.join(random.choices(string.hexdigits, k=8))}


class TestPublishRevision(BaseRevisionTest):
    """Test case for publishing a revision."""

    def setUp(self):
        super().setUp()
        self.base_path = "/dp/v1/collections"
        self.mock_timestamp = datetime(2000, 12, 25, 0, 0)


    @patch("backend.corpora.common.entities.dataset_asset.s3_client.head_object", wraps=get_random_etag)
    def publish_collection(self, collection_id: str, mocked_func) -> dict:
        """
        Verify publish a collection under revision.
        :return: Jsonified response of GET collection/<collection_id>.
        """
        # Check initial published_at and revised_at. Since collection creation
        # for the already published collection/datasets do not go through the
        # normal user flow, no initial values for published_at and revised_at
        self.assertIsNone(self.pub_collection.published_at)
        self.assertIsNone(self.pub_collection.revised_at)

        for dataset in self.pub_collection.datasets:
            self.assertIsNone(dataset.published_at)
            self.assertIsNone(dataset.revised_at)

        self.session.expire_all()

        path = f"{self.base_path}/{collection_id}/publish"
        with patch("backend.corpora.common.entities.collection.datetime") as mock_dt:
            mock_dt.utcnow = Mock(return_value=self.mock_timestamp)
            response = self.app.post(path, headers=self.headers)
        self.assertEqual(202, response.status_code)

        self.assertDictEqual({"collection_uuid": collection_id, "visibility": "PUBLIC"}, json.loads(response.data))
        self.addCleanup(self.delete_collection, collection_id, "PUBLIC")

        # Cannot call publish for an already published collection
        response = self.app.post(path, headers=self.headers)
        self.assertEqual(403, response.status_code)

        # Check that the published collection is listed in /collections
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(self.base_path, headers=headers)
        self.assertEqual(200, response.status_code)

        ids = [col["id"] for col in json.loads(response.data)["collections"]]
        self.assertIn(collection_id, ids)

        # Check GET collection/<collection_id>
        path = f"{self.base_path}/{collection_id}"
        response = self.app.get(path, headers=headers)
        self.assertEqual(200, response.status_code)

        response_json = json.loads(response.data)
        self.assertEqual("PUBLIC", response_json["visibility"])
        self.assertEqual(collection_id, response_json["id"])

        return response_json

    def verify_datasets(self, actual_body: dict, expected_dataset_ids: typing.List[str]) -> None:
        """Verify collection datasets."""
        actual_datasets = [d["id"] for d in actual_body["datasets"]]
        self.assertListEqual(expected_dataset_ids, actual_datasets)
        self.assertTrue(all([d["published"] for d in actual_body["datasets"]]))

        for dataset_id in expected_dataset_ids:
            dataset = Dataset.get(self.session, dataset_id)
            self.assertIn(dataset.id, dataset.explorer_url)
            for s3_object in self.get_s3_object_paths_from_dataset(dataset):
                if dataset.tombstone:
                    self.assertS3FileDoesNotExist(*s3_object)
                else:
                    self.assertS3FileExists(*s3_object)

    def test__with_revision_with_new_dataset__OK(self):
        """Publish a revision with new datasets."""
        new_dataset_id = self.generate_dataset_with_s3_resources(
            self.session, collection_id=self.rev_collection.id, collection_visibility=CollectionVisibility.PRIVATE
        ).id
        dataset_ids = [ds.id for ds in self.pub_collection.datasets]
        dataset_ids.append(new_dataset_id)

        # Publish revision
        response_json = self.publish_collection(self.rev_collection.id)
        self.verify_datasets(response_json, dataset_ids)

        # Check published_at and revised_at
        # Collection: Only revised_at should be updated
        self.assertIsNone(response_json.get("published_at"))
        self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(response_json["revised_at"]))

        # Datasets: Only the newly added dataset should have published_at updated
        for dataset in response_json["datasets"]:
            if dataset["id"] == new_dataset_id:
                self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(dataset["published_at"]))
            else:
                self.assertIsNone(dataset.get("published_at"))
            self.assertIsNone(dataset.get("revised_at"))

    def test__with_revision_with_tombstoned_datasets__OK(self):
        """Publish a revision with delete datasets."""
        rev_dataset_id = self.rev_collection.datasets[0].id
        pub_dataset = self.pub_collection.datasets[0]
        published_s3_objects = self.get_s3_object_paths_from_dataset(pub_dataset)

        # Delete a dataset under revision
        self.app.delete(f"/dp/v1/datasets/{rev_dataset_id}", headers=self.headers)

        # Publish the revision with the deleted dataset
        response_json = self.publish_collection(self.rev_collection.id)
        self.session.expire_all()

        dataset = Dataset.get(self.session, pub_dataset.id, include_tombstones=True)

        self.assertIn(pub_dataset.id, dataset.explorer_url)
        self.assertTrue(dataset.tombstone)
        for s3_object in published_s3_objects:
            self.assertS3FileDoesNotExist(*s3_object)

        # Check published_at and revised_at
        # Collection: Only revised_at should be updated
        self.assertIsNone(response_json.get("published_at"))
        self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(response_json["revised_at"]))

        # Datasets: None should be updated
        self.assertIsNone(dataset.published_at)
        self.assertIsNone(dataset.revised_at)
        for dataset in response_json["datasets"]:
            self.assertIsNone(dataset.get("published_at"))
            self.assertIsNone(dataset.get("revised_at"))

    def test__with_revision_with_tombstoned_datasets_rollback__OK(self):
        """
        Revision state is restored and s3 assets are unchanged, if the database
        transactions fails.
        """
        rev_dataset_id = self.rev_collection.datasets[0].id
        pub_dataset = self.pub_collection.datasets[0]
        published_s3_objects = self.get_s3_object_paths_from_dataset(pub_dataset)

        # Delete a dataset under revision
        self.app.delete(f"/dp/v1/datasets/{rev_dataset_id}", self.headers)
        self.session.expire_all()

        # Publish revision with database transaction failure
        expect_collection = self.rev_collection.reshape_for_api(tombstoned_datasets=True)
        expect_datasets = expect_collection.pop("datasets")
        expect_datasets.sort(key=lambda x: x["id"])
        with self.assertRaises(CorporaException):
            with db_session_manager() as session:
                rev_collection = Collection.get(session, (self.rev_collection.id, self.rev_collection.visibility))
                with mock.patch.object(rev_collection.session, "commit", side_effect=SQLAlchemyError):
                    rev_collection.publish()

        self.session.expire_all()

        actual_collection = self.rev_collection.reshape_for_api(tombstoned_datasets=True)
        actual_datasets = actual_collection.pop("datasets")
        actual_datasets.sort(key=lambda x: x["id"])
        self.assertEqual(expect_collection, actual_collection)
        self.assertEqual(expect_datasets, actual_datasets)

        dataset = Dataset.get(self.session, pub_dataset.id, include_tombstones=True)
        self.assertFalse(dataset.tombstone)
        for s3_object in published_s3_objects:
            self.assertS3FileExists(*s3_object)

        # Check published_at and revised_at
        # Collection: None should be updated
        self.assertIsNone(actual_collection.get("published_at"))
        self.assertIsNone(actual_collection.get("revised_at"))

        # Datasets: None should be updated
        self.assertIsNone(dataset.published_at)
        self.assertIsNone(dataset.revised_at)
        for dataset in actual_datasets:
            self.assertIsNone(dataset.get("published_at"))
            self.assertIsNone(dataset.get("revised_at"))

    def test__with_revision_with_all_tombstoned_datasets__409(self):
        """Unable to publish a revision with no datasets."""
        for dataset in self.rev_collection.datasets:
            self.app.delete(f"/dp/v1/datasets/{dataset.id}", headers=self.headers)
        path = f"/dp/v1/collections/{self.rev_collection.id}/publish"
        response = self.app.post(path, headers=self.headers)
        self.assertEqual(409, response.status_code)

    def test__with_revision_with_refreshed_datasets__OK(self):
        """"Publish a revision with refreshed datasets."""
        self.refresh_datasets()

        response_json = self.publish_collection(self.rev_collection.id)
        self.verify_datasets(response_json, [ds.id for ds in self.pub_collection.datasets])

        # Check published_at and revised_at
        # Collection: Only revised_at should be updated
        self.assertIsNone(response_json.get("published_at"))
        self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(response_json["revised_at"]))

        # Datatsets: None should be updated
        for dataset in response_json["datasets"]:
            self.assertIsNone(dataset.get("published_at"))
            self.assertEqual(self.mock_timestamp, datetime.utcfromtimestamp(dataset["revised_at"]))

    def test__publish_revision_with_collection_info_updated__201(self):
        """Publish a revision with collection detail changes."""
        expected_body = {
            "name": "collection name",
            "description": "This is a test collection",
            "contact_name": "person human",
            "contact_email": "person@human.com",
            "links": [
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                {"link_name": "DOI Link 2", "link_url": "http://doi.org/10.1017", "link_type": "DOI"},
            ],
            "data_submission_policy_version": "0",
        }
        self.rev_collection.update(**expected_body)
        pub_s3_objects, _ = self.get_s3_objects_from_collections()

        # Published revision with collection details updated
        response_json = self.publish_collection(self.rev_collection.id)
        self.assertPublishedCollectionOK(expected_body, pub_s3_objects)

        self.verify_datasets(response_json, [ds.id for ds in self.pub_collection.datasets])

        # Check published_at and revised_at
        # Collection: None should be updated
        self.assertIsNone(response_json.get("published_at"))
        self.assertIsNone(response_json.get("revised_at"))

        # Datasets: None should be updated
        for dataset in response_json["datasets"]:
            self.assertIsNone(dataset.get("published_at"))
            self.assertIsNone(dataset.get("revised_at"))

    def test__with_revision_and_existing_datasets(self):
        """Publish a revision with the same, existing datasets."""
        response_json = self.publish_collection(self.rev_collection.id)
        self.verify_datasets(response_json, [ds.id for ds in self.pub_collection.datasets])

        # Check published_at and revised_at
        # Collection: None should be updated
        self.assertIsNone(response_json.get("published_at"))
        self.assertIsNone(response_json.get("revised_at"))

        # Datasets: None should be updated
        for dataset in response_json["datasets"]:
            self.assertIsNone(dataset.get("published_at"))
            self.assertIsNone(dataset.get("revised_at"))
