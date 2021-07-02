import json
import typing
from unittest import mock

from sqlalchemy.exc import SQLAlchemyError

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.entities import Dataset, Collection
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.exceptions import CorporaException
from backend.corpora.common.utils.json import CustomJSONEncoder
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token


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
            self.assertEqual(resp.status_code, 200)
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
    def test__start_revision__201(self):
        tests = [
            {"visibility": CollectionVisibility.PUBLIC},
            {
                "visibility": CollectionVisibility.PUBLIC,
                "links": [
                    {"link_name": "Link 1", "link_url": "This is a new link", "link_type": "OTHER"},
                    {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
                ],
            },
        ]
        for test in tests:
            with self.subTest(test):
                collection = self.generate_collection(self.session, **test)
                dataset_0 = self.generate_dataset_with_s3_resources(self.session, collection=collection, published=True)
                dataset_1 = self.generate_dataset_with_s3_resources(self.session, collection=collection, published=True)

                test_url = f"/dp/v1/collections/{collection.id}"
                headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

                # Test post
                response = self.app.post(test_url, headers=headers)
                self.assertEqual(201, response.status_code)
                post_body = json.loads(response.data)
                for key in test.keys():
                    if key == "visibility":
                        self.assertEqual("PRIVATE", post_body[key])
                    else:
                        self.assertEqual(test[key], post_body[key])
                # Test get
                response = self.app.get(f"{test_url}?visibility=PRIVATE", headers=headers)
                self.assertEqual(200, response.status_code)
                get_body = json.loads(response.data)
                self.assertEqual(post_body, get_body)
                with self.subTest("Test datasets in revised collection are not original datasets"):
                    new_dataset_ids = [x["id"] for x in get_body["datasets"]]
                    self.assertNotIn(dataset_0.id, new_dataset_ids)
                    self.assertNotIn(dataset_1.id, new_dataset_ids)
                with self.subTest("Test revised datasets point at original datasets"):
                    original_dataset_ids = [x["original_id"] for x in get_body["datasets"]]
                    self.assertIn(dataset_0.id, original_dataset_ids)
                    self.assertIn(dataset_1.id, original_dataset_ids)
                with self.subTest("Check assets point at revised dataset"):
                    assets_0 = get_body["datasets"][0]["dataset_assets"]
                    assets_1 = get_body["datasets"][1]["dataset_assets"]
                    revised_dataset_0 = get_body["datasets"][0]["id"]
                    revised_dataset_1 = get_body["datasets"][1]["id"]
                    for x in assets_0:
                        self.assertEqual(revised_dataset_0, x["dataset_id"])
                    for x in assets_1:
                        self.assertEqual(revised_dataset_1, x["dataset_id"])
                # Test unauthenticated get
                get_body.pop("access_type")
                expected_body = get_body
                headers = {"host": "localhost", "Content-Type": "application/json"}
                response = self.app.get(f"{test_url}?visibility=PRIVATE", headers=headers)
                self.assertEqual(200, response.status_code)
                actual_body = json.loads(response.data)
                self.assertEqual("READ", actual_body.pop("access_type"))
                self.assertEqual(expected_body, actual_body)

    def test__revision__409(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC)
        test_url = f"/dp/v1/collections/{collection.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(201, response.status_code)
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(409, response.status_code)

    def test__revision_nonexistent__403(self):
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.post("/dp/v1/collections/random", headers=headers)
        self.assertEqual(403, response.status_code)

    def test__revision_not_owner__403(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PUBLIC, owner="someone else"
        )
        test_url = f"/dp/v1/collections/{collection.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        response = self.app.post(test_url, headers=headers)
        self.assertEqual(403, response.status_code)


class TestDeleteRevision(BaseRevisionTest):
    def setUp(self):
        super().setUp()
        url = f"/dp/v1/collections/{self.pub_collection.id}"
        self.test_url_collect_private = f"{url}?visibility=PRIVATE"
        self.test_url_collection_public = f"{url}?visibility=PUBLIC"

    def test__revision_deleted__204(self):
        # Delete the revision
        resp = self.app.delete(self.test_url_collect_private, headers=self.headers)
        self.assertEqual(204, resp.status_code)

        # Cannot get the revision
        resp = self.app.get(self.test_url_collect_private, headers=self.headers)
        self.assertEqual(403, resp.status_code)

    def test__revision_deleted_with_published_datasets(self):
        """The published dataset artifacts should be intact after deleting a collection revision"""
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
        """The refreshed datasets should be deleted and the published dataset intact. The published dataset artifacts should
        be intact after deleting a collection revision
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
        """The refreshed datasets should be deleted and the published dataset restored in the revision."""
        expected_pub_body = json.loads(json.dumps(self.pub_collection.reshape_for_api(), cls=CustomJSONEncoder))
        expected_rev_body = self.remove_timestamps(
            json.loads(json.dumps(self.rev_collection.reshape_for_api(), cls=CustomJSONEncoder)),
            remove=["id", "dataset_id"],
        )
        refreshed_datasets = self.refresh_datasets()
        pub_s3_objects, rev_s3_objects = self.get_s3_objects_from_collections()

        # Refreshed datasets do not point to the published resources in s3.
        for s3_object in rev_s3_objects:
            self.assertNotIn(s3_object, pub_s3_objects)

        # Delete the refreshed_datasets
        for ds in refreshed_datasets:
            resp = self.app.delete(f"/dp/v1/datasets/{ds.id}", headers=self.headers)
            self.assertEqual(resp.status_code, 202)

        with self.subTest("refreshed artifacts and explorer s3 objects deleted"):
            for bucket, file_name in rev_s3_objects:
                self.assertS3FileDoesNotExist(bucket, file_name)
        self.assertPublishedCollectionOK(expected_pub_body, pub_s3_objects)

        # Check that the original datasets info has been restored to the revision dataset.
        self.session.expire_all()
        actual_rev_body = self.remove_timestamps(
            json.loads(json.dumps(self.rev_collection.reshape_for_api(), cls=CustomJSONEncoder)),
            remove=["id", "dataset_id"],
        )
        self.assertEqual(expected_rev_body, actual_rev_body)

    def test__delete_published_dataset_during_revision(self):
        """The dataset is tombstone in the revision. The published artifacts are intact"""
        expected_body = json.loads(json.dumps(self.pub_collection.reshape_for_api(), cls=CustomJSONEncoder))
        pub_s3_objects, rev_s3_objects = self.get_s3_objects_from_collections()
        # Delete a published dataset in the revision
        rev_dataset_count = len(self.rev_collection.datasets)
        rev_dataset_id = self.rev_collection.datasets[0].id
        test_dataset_url = f"/dp/v1/datasets/{rev_dataset_id}"
        resp = self.app.delete(test_dataset_url, headers=self.headers)
        self.assertEqual(resp.status_code, 202)

        # Get the revision authenticated
        resp = self.app.get(self.test_url_collect_private, headers=self.headers)
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(rev_dataset_count, len(json.loads(resp.data)["datasets"]))

        # Get the revision unauthenticated
        resp = self.app.get(self.test_url_collect_private)
        self.assertEqual(resp.status_code, 200)
        # The dataset is a tombstone in the revision
        self.assertEqual(rev_dataset_count - 1, len(json.loads(resp.data)["datasets"]))
        self.session.expire_all()
        for dataset in self.rev_collection.datasets:
            if dataset.id == rev_dataset_id:
                self.assertTrue(dataset.tombstone)
                break
        self.assertPublishedCollectionOK(expected_body, pub_s3_objects)


class TestPublishRevision(BaseRevisionTest):
    def setUp(self):
        super().setUp()
        self.test_url = f"/dp/v1/collections/{self.pub_collection.id}/publish"

    def publish_collection(self, collection_id: str):
        self.session.expire_all()
        path = f"/dp/v1/collections/{collection_id}/publish"
        response = self.app.post(path, headers=self.headers)
        self.assertEqual(202, response.status_code)
        self.assertDictEqual({"collection_uuid": collection_id, "visibility": "PUBLIC"}, json.loads(response.data))
        self.addCleanup(self.delete_collection, collection_id, "PUBLIC")

        # cannot call twice
        response = self.app.post(path, headers=self.headers)
        self.assertEqual(403, response.status_code)

        # check if the collection is listed
        path = "/dp/v1/collections"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(path, headers)
        self.assertEqual(200, response.status_code)
        ids = [col["id"] for col in json.loads(response.data)["collections"]]
        self.assertIn(collection_id, ids)

        # check get collection_uuid
        path = f"/dp/v1/collections/{collection_id}"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        response = self.app.get(path, headers)
        self.assertEqual(200, response.status_code)
        actual = json.loads(response.data)
        self.assertEqual("PUBLIC", actual["visibility"])
        self.assertEqual(collection_id, actual["id"])
        return actual

    def verify_datasets(self, actual_body, expected_dataset_ids):
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
        """publish a revision with new datasets"""
        new_dataset_id = self.generate_dataset_with_s3_resources(
            self.session, collection_id=self.rev_collection.id, collection_visibility=CollectionVisibility.PRIVATE
        ).id
        dataset_ids = [ds.id for ds in self.pub_collection.datasets]
        dataset_ids.append(new_dataset_id)
        actual_body = self.publish_collection(self.rev_collection.id)
        self.verify_datasets(actual_body, [ds.id for ds in self.pub_collection.datasets])

    def test__with_revision_with_tombstoned_datasets__OK(self):
        """publish a revision with delete datasets"""
        rev_dataset_id = self.rev_collection.datasets[0].id
        pub_dataset = self.pub_collection.datasets[0]
        published_s3_objects = self.get_s3_object_paths_from_dataset(pub_dataset)
        self.app.delete(f"/dp/v1/datasets/{rev_dataset_id}", headers=self.headers)
        self.publish_collection(self.rev_collection.id)
        self.session.expire_all()
        dataset = Dataset.get(self.session, pub_dataset.id, include_tombstones=True)
        self.assertIn(pub_dataset.id, dataset.explorer_url)
        self.assertTrue(dataset.tombstone)
        for s3_object in published_s3_objects:
            self.assertS3FileDoesNotExist(*s3_object)

    def test__with_revision_with_tombstoned_datasets_rollback__OK(self):
        """revision state is restored and s3 assets are unchanged, if the database transactions fails."""
        rev_dataset_id = self.rev_collection.datasets[0].id
        pub_dataset = self.pub_collection.datasets[0]
        published_s3_objects = self.get_s3_object_paths_from_dataset(pub_dataset)

        self.app.delete(f"/dp/v1/datasets/{rev_dataset_id}", self.headers)
        self.session.expire_all()

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

    def test__with_revision_with_all_tombstoned_datasets__409(self):
        """unable to publish a revision with no datasets"""
        for dataset in self.rev_collection.datasets:
            self.app.delete(f"/dp/v1/datasets/{dataset.id}", headers=self.headers)
        path = f"/dp/v1/collections/{self.rev_collection.id}/publish"
        response = self.app.post(path, headers=self.headers)
        self.assertEqual(409, response.status_code)

    def test__with_revision_with_refreshed_datasets__OK(self):
        """"publish a revision with refreshed datasets"""
        self.refresh_datasets()
        actual_body = self.publish_collection(self.rev_collection.id)
        self.verify_datasets(actual_body, [ds.id for ds in self.pub_collection.datasets])

    def test__publish_revision_with_collection_info_updated__201(self):
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
        pub_s3_objects, rev_s3_objects = self.get_s3_objects_from_collections()
        actual_body = self.publish_collection(self.rev_collection.id)
        self.assertPublishedCollectionOK(expected_body, pub_s3_objects)
        self.verify_datasets(actual_body, [ds.id for ds in self.pub_collection.datasets])

    def test_with_revision_and_existing_datasets(self):
        actual_body = self.publish_collection(self.rev_collection.id)
        self.verify_datasets(actual_body, [ds.id for ds in self.pub_collection.datasets])
