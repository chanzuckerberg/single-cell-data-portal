import botocore
import json
from backend.corpora.common.corpora_orm import CollectionVisibility, DatasetArtifactFileType
from backend.corpora.common.entities import DatasetAsset
from backend.corpora.common.entities.dataset import get_cxg_bucket_path
from backend.corpora.common.utils.db_session import clone
from backend.corpora.common.utils.json import CustomJSONEncoder
from tests.unit.backend.chalice.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from tests.unit.backend.chalice.api_server.mock_auth import get_auth_token


class TestRevision(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def test__revision__201(self):
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

    def test__publish_revision__201(self):
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

    def test__revision_deleted_with_published_datasets(self):
        """The published dataset artifacts should be intact after deleting a collection revision"""
        # Generate public collection
        pub_collection = self.generate_collection(self.session, visibility="PUBLIC")
        pub_dataset = self.generate_dataset(
            self.session, collection_visibility="PUBLIC", collection_id=pub_collection.id, published=True
        )
        dep_dir = self.generate_deployment_directory(self.session, pub_dataset.id, upload=True)
        arts = [
            self.generate_artifact(self.session, pub_dataset.id, ext, upload=True) for ext in DatasetArtifactFileType
        ]
        expected_body = json.loads(json.dumps(pub_collection.reshape_for_api(), cls=CustomJSONEncoder))

        # Generate revision
        rev_collection = self.generate_collection(self.session, visibility="PRIVATE", id=pub_collection.id)
        rev_datasets = self.generate_dataset(
            self.session,
            collection_visibility="PRIVATE",
            collection_id=rev_collection.id,
            published=True,
            original_id=pub_dataset.id,
        )
        for child in arts:
            self.session.add(clone(child.db_object, dataset_id=rev_datasets.id))
        self.session.add(clone(dep_dir, dataset_id=rev_datasets.id))
        self.session.commit()

        test_url = f"/dp/v1/collections/{pub_collection.id}?visibility=PRIVATE"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        resp = self.app.delete(test_url, headers=headers)
        resp.raise_for_status()
        resp = self.app.get(test_url, headers=headers)
        self.assertEqual(403, resp.status_code)

        with self.subTest("published artifacts ok"):
            for artifact in pub_dataset.artifacts:
                art = DatasetAsset(artifact)
                self.assertGreater(art.get_file_size(), 1)
        with self.subTest("published deployed directories ok"):
            s3_file = get_cxg_bucket_path(pub_dataset.deployment_directories[0])
            self.assertGreater(self.cellxgene_bucket.Object(s3_file).content_length, 1)
        with self.subTest("publish collection ok"):
            test_url = f"/dp/v1/collections/{pub_collection.id}"
            resp = self.app.get(test_url)
            resp.raise_for_status()
            actual_body = json.loads(resp.body)
            for key in expected_body.keys():
                if key in actual_body.keys():
                    self.assertEqual(expected_body[key], actual_body[key])

    def test__revision_deleted_with_new_datasets(self):
        """The new datasets should be deleted when the revison is deleted."""
        # Generate public collection
        pub_collection = self.generate_collection(self.session, visibility="PUBLIC")

        # Generate revision
        rev_collection = self.generate_collection(self.session, visibility="PRIVATE", id=pub_collection.id)
        rev_datasets = self.generate_dataset(
            self.session,
            collection_visibility="PRIVATE",
            collection_id=rev_collection.id,
            published=False,
        )
        self.generate_deployment_directory(self.session, rev_datasets.id, upload=True)
        [self.generate_artifact(self.session, rev_datasets.id, ext, upload=True) for ext in DatasetArtifactFileType]

        with self.subTest("new artifacts exist"):
            for artifact in rev_datasets.artifacts:
                art = DatasetAsset(artifact)
                self.assertGreater(art.get_file_size(), 1)
        with self.subTest("new deployed directories exist"):
            s3_file = get_cxg_bucket_path(rev_datasets.deployment_directories[0])
            self.assertGreater(self.cellxgene_bucket.Object(s3_file).content_length, 1)

        test_url = f"/dp/v1/collections/{pub_collection.id}?visibility=PRIVATE"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        resp = self.app.delete(test_url, headers=headers)
        resp.raise_for_status()
        resp = self.app.get(test_url, headers=headers)
        self.assertEqual(403, resp.status_code)

        with self.subTest("new artifacts deleted"):
            for artifact in rev_datasets.artifacts:
                art = DatasetAsset(artifact)
                self.assertIsNone(art.get_file_size())
        with self.subTest("new deployed directories deleted"):
            with self.assertRaises(botocore.exceptions.ClientError):
                s3_file = get_cxg_bucket_path(rev_datasets.deployment_directories[0])
                self.assertGreater(self.cellxgene_bucket.Object(s3_file).content_length, 1)

    def test__revision_deleted_with_refreshed_datasets(self):
        """The refreshed datasets should be deleted and the published dataset intact."""
        """The published dataset artifacts should be intact after deleting a collection revision"""
        # Generate public collection
        pub_collection = self.generate_collection(self.session, visibility="PUBLIC")
        pub_dataset = self.generate_dataset(
            self.session, collection_visibility="PUBLIC", collection_id=pub_collection.id, published=True
        )
        self.generate_deployment_directory(self.session, pub_dataset.id, upload=True)
        [self.generate_artifact(self.session, pub_dataset.id, ext, upload=True) for ext in DatasetArtifactFileType]
        expected_body = json.loads(json.dumps(pub_collection.reshape_for_api(), cls=CustomJSONEncoder))

        # Generate revision
        rev_collection = self.generate_collection(self.session, visibility="PRIVATE", id=pub_collection.id)
        rev_datasets = self.generate_dataset(
            self.session,
            collection_visibility="PRIVATE",
            collection_id=rev_collection.id,
            published=False,
            original_id=pub_dataset.id,
        )
        dep_dir = self.generate_deployment_directory(self.session, rev_datasets.id, upload=True)
        arts = [
            self.generate_artifact(self.session, rev_datasets.id, ext, upload=True) for ext in DatasetArtifactFileType
        ]
        s3_objects = [(self.bucket, art.key_name) for art in arts] + [
            (self.cellxgene_bucket, get_cxg_bucket_path(dep_dir))
        ]

        test_url = f"/dp/v1/collections/{pub_collection.id}?visibility=PRIVATE"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        resp = self.app.delete(test_url, headers=headers)
        resp.raise_for_status()
        resp = self.app.get(test_url, headers=headers)
        self.assertEqual(403, resp.status_code)

        with self.subTest("published artifacts ok"):
            for artifact in pub_dataset.artifacts:
                art = DatasetAsset(artifact)
                self.assertGreater(art.get_file_size(), 1)
        with self.subTest("published deployed directories ok"):
            s3_file = get_cxg_bucket_path(pub_dataset.deployment_directories[0])
            self.assertGreater(self.cellxgene_bucket.Object(s3_file).content_length, 1)
        with self.subTest("publish collection ok"):
            test_url = f"/dp/v1/collections/{pub_collection.id}"
            resp = self.app.get(test_url)
            resp.raise_for_status()
            actual_body = json.loads(resp.body)
            for key in expected_body.keys():
                if key in actual_body.keys():
                    self.assertEqual(expected_body[key], actual_body[key])
        with self.subTest("refreshed dataset s3 resources deleted"):
            for bucket, key in s3_objects:
                with self.assertRaises(botocore.exceptions.ClientError):
                    bucket.Object(key).content_length
