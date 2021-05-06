import json
import typing

from backend.corpora.common.corpora_orm import CollectionVisibility
from backend.corpora.common.utils.db_session import clone
from backend.corpora.common.utils.json import CustomJSONEncoder
from tests.unit.backend.chalice.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from tests.unit.backend.chalice.api_server.mock_auth import get_auth_token


class TestRevision(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
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
                test_url = f"/dp/v1/collections/{collection.id}"
                headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

                # Test post
                response = self.app.post(test_url, headers=headers)
                self.assertEqual(201, response.status_code)
                post_body = json.loads(response.body)
                for key in test.keys():
                    if key == "visibility":
                        self.assertEqual("PRIVATE", post_body[key])
                    else:
                        self.assertEqual(test[key], post_body[key])
                # Test get
                response = self.app.get(f"{test_url}?visibility=PRIVATE", headers=headers)
                self.assertEqual(200, response.status_code)
                get_body = json.loads(response.body)
                self.assertEqual(post_body, get_body)

                # Test unauthenticated get
                get_body.pop("access_type")
                expected_body = get_body
                headers = {"host": "localhost", "Content-Type": "application/json"}
                response = self.app.get(f"{test_url}?visibility=PRIVATE", headers=headers)
                self.assertEqual(200, response.status_code)
                actual_body = json.loads(response.body)
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

    def test__publish_revision_with_collection_info_updated__201(self):
        collection = self.generate_collection(
            self.session,
            visibility=CollectionVisibility.PUBLIC,
            links=[
                {"link_url": "http://doi.org/10.1010", "link_type": "DOI"},
                {"link_name": "DOI Link", "link_url": "http://doi.org/10.1016", "link_type": "DOI"},
            ],
        )
        test_url = f"/dp/v1/collections/{collection.id}"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

        # Start revision
        response = self.app.post(test_url, headers=headers)
        response.raise_for_status()

        # Update the revision
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
        data = json.dumps(expected_body)
        response = self.app.put(f"{test_url}?visibility=PRIVATE", data=data, headers=headers)
        response.raise_for_status()

        # publish the revision
        response = self.app.post(f"{test_url}/publish", headers=headers)
        response.raise_for_status()

        # Get the revised collection
        response = self.app.get(f"{test_url}", headers=headers)
        response.raise_for_status()
        actual_body = json.loads(response.body)
        self.assertEqual(actual_body["visibility"], "PUBLIC")
        for link in expected_body.pop("links"):
            self.assertIn(link, actual_body["links"])
        for keys in expected_body.keys():
            self.assertEqual(expected_body[keys], actual_body[keys])


class TestDeleteRevision(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        super().setUp()
        pub_collection = self.generate_collection(self.session, visibility="PUBLIC")
        self.generate_dataset_with_s3_resources(
            self.session, collection_visibility="PUBLIC", collection_id=pub_collection.id, published=True
        )
        url = f"/dp/v1/collections/{pub_collection.id}"
        self.test_url_collect_private = f"{url}?visibility=PRIVATE"
        self.test_url_collection_public = f"{url}?visibility=PUBLIC"
        self.pub_collection = pub_collection
        self.rev_collection = pub_collection.revision()
        self.headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

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
        pub_s3_objects, rev_s3_objects = self.copy_dataset_from_published_to_revision()

        # Revision and Published collection refer to the same S3 resources
        self.assertEqual(pub_s3_objects, rev_s3_objects)

        # Delete the revision
        resp = self.app.delete(self.test_url_collect_private, headers=self.headers)
        resp.raise_for_status()
        self.assertPublishedCollectionOK(expected_body, pub_s3_objects)

    def test__revision_deleted_with_new_datasets(self):
        """The new datasets should be deleted when the revison is deleted."""
        # Generate revision dataset
        rev_dataset = self.generate_dataset_with_s3_resources(
            self.session,
            collection_visibility="PRIVATE",
            collection_id=self.rev_collection.id,
            published=False,
        )
        s3_objects = self.get_s3_object_paths_from_dataset(rev_dataset)

        # Check resources exist
        with self.subTest("new artifacts and deployment_directories exist"):
            for bucket, file_name in s3_objects:
                self.assertS3FileExists(bucket, file_name)

        # Delete Revision
        resp = self.app.delete(self.test_url_collect_private, headers=self.headers)
        resp.raise_for_status()

        with self.subTest("new artifacts and deployment_directories deleted"):
            for bucket, file_name in s3_objects:
                self.assertS3FileDoesNotExist(bucket, file_name)

    def test__revision_deleted_with_refreshed_datasets(self):
        """The refreshed datasets should be deleted and the published dataset intact. The published dataset artifacts should
        be intact after deleting a collection revision
        """
        expected_body = json.loads(json.dumps(self.pub_collection.reshape_for_api(), cls=CustomJSONEncoder))
        pub_s3_objects, rev_s3_objects = self.copy_dataset_from_published_to_revision(refreshed=True)

        # Refreshed datasets do not point to the published resources in s3.
        for s3_object in rev_s3_objects:
            self.assertNotIn(s3_object, pub_s3_objects)

        # Delete the revision
        resp = self.app.delete(self.test_url_collect_private, headers=self.headers)
        resp.raise_for_status()

        with self.subTest("refreshed artifacts and deployment_directories deleted"):
            for bucket, file_name in rev_s3_objects:
                self.assertS3FileDoesNotExist(bucket, file_name)
        self.assertPublishedCollectionOK(expected_body, pub_s3_objects)

    def test__delete_published_dataset_during_revision(self):
        """The dataset is tombstone in the revision. The published artifacts are intact"""
        # Generate public collection
        expected_body = json.loads(json.dumps(self.pub_collection.reshape_for_api(), cls=CustomJSONEncoder))
        pub_s3_objects, rev_s3_objects = self.copy_dataset_from_published_to_revision()
        # Delete a published dataset in the revision
        rev_dataset_id = self.rev_collection.datasets[0].id
        test_dataset_url = f"/dp/v1/datasets/{rev_dataset_id}"
        resp = self.app.delete(test_dataset_url, headers=self.headers)
        resp.raise_for_status()

        # Get the revision
        resp = self.app.get(self.test_url_collect_private, headers=self.headers)
        resp.raise_for_status()

        # The dataset is a tombstone in the revision
        self.assertEqual(json.loads(resp.body)["datasets"], [])
        self.session.expire_all()
        for dataset in self.rev_collection.datasets:
            if dataset.id == rev_dataset_id:
                self.assertTrue(self.rev_collection.datasets[0].tombstone)
                break
        self.assertPublishedCollectionOK(expected_body, pub_s3_objects)

    def assertPublishedCollectionOK(self, expected_body, s3_objects):
        with self.subTest("published artifacts and deployment_directories ok"):
            for bucket, file_name in s3_objects:
                self.assertS3FileExists(bucket, file_name)
        with self.subTest("publish collection ok"):
            resp = self.app.get(self.test_url_collection_public)
            resp.raise_for_status()
            actual_body = json.loads(resp.body)
            for key in expected_body.keys():
                if key in actual_body.keys():
                    self.assertEqual(expected_body[key], actual_body[key])

    def copy_dataset_from_published_to_revision(self, refreshed=False) -> typing.Tuple[typing.List, typing.List]:
        """

        :param refreshed: created a refreshed version of the published dataset
        :return: a list of s3 objects in the published collection, and a list of s3 objects the revision collection.
        """
        rev_s3_objects = []
        pub_s3_objects = []

        # Copy published datasets to revision
        for dataset in self.pub_collection.datasets:
            if refreshed:
                rev_dataset = self.generate_dataset_with_s3_resources(
                    self.session,
                    collection_visibility="PRIVATE",
                    collection_id=self.rev_collection.id,
                    original_id=dataset.id,
                )
            else:
                rev_dataset = self.generate_dataset(
                    self.session,
                    collection_visibility="PRIVATE",
                    collection_id=self.rev_collection.id,
                    published=True,
                    original_id=dataset.id,
                )
                for artifact in dataset.artifacts:
                    self.session.add(clone(artifact, dataset_id=rev_dataset.id))
                self.session.add(clone(dataset.deployment_directories[0], dataset_id=rev_dataset.id))
                self.session.commit()
            pub_s3_objects.extend(self.get_s3_object_paths_from_dataset(dataset))
            rev_s3_objects.extend(self.get_s3_object_paths_from_dataset(rev_dataset))
        return pub_s3_objects, rev_s3_objects
