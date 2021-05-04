import json
from backend.corpora.common.corpora_orm import CollectionVisibility, DatasetArtifactFileType, DatasetArtifactType
from backend.corpora.common.utils.db_session import clone
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

    # test a revision can be delete
    # test with published datasets. The published datasets should not be delete if the revision is canceled
    # def test__revision__deleted(self):
    #     """ with a published datasets"""
    #     # Generate public collection
    #     pub_collection = self.generate_collection(self.session, visibility="PUBLIC")
    #     pub_datasets = self.generate_dataset(
    #         self.session, collection_visibility="PUBLIC", collection_id=pub_collection.id, published=True
    #     )
    #     for ds in pub_datasets:
    #         # create artifacts and upload to s3
    #         ds.update(
    #             artifacts=[self.create_artifact(ds.id, ext.name) for ext in DatasetArtifactFileType],
    #             deployment_directories=[dict(url=f"http://bogus.url/d/{ds.id}.cxg/")],
    #         )
    #     # Generate revision
    #     rev_collection = self.generate_collection(self.session, visibility="PRIVATE", id=pub_collection.id)
    #     rev_datasets = [
    #         self.generate_dataset(
    #             self.session,
    #             collection_visibility="PRIVATE",
    #             collection_id=rev_collection.id,
    #             published=True,
    #             original_id=pub_datasets[0].id,
    #         ),
    #         self.generate_dataset(
    #             self.session,
    #             collection_visibility="PRIVATE",
    #             collection_id=rev_collection.id,
    #             published=False,
    #             original_id=pub_datasets[1].id,
    #         ),
    #         self.generate_dataset(
    #             self.session, collection_visibility="PRIVATE", collection_id=rev_collection.id, published=False
    #         ),
    #     ]
    #     for ds in rev_datasets:
    #         # Clone the artifacts
    #         ds.update(
    #             artifacts=[self.create_artifact(ds.id, ext.name, False) for ext in DatasetArtifactFileType],
    #             deployment_directories=[dict(url=f"http://bogus.url/d/{ds.id}.cxg/")],
    #         )
    #     test_url = f"/dp/v1/collections/{pub_collection.id}?visibility=PRIVATE"
    #     headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
    #     self.app.delete(test_url, headers=headers)
    #     for dataset in pub_datasets:
    #         for artifact in dataset.artifacts:
    #             with self.subTest(artifact.filename):
    #                 self.assertGreater(self.bucket.Object(artifact.filename).content_length, 1)

    def test__revision__deleted(self):
        """ with a published datasets"""
        # Generate public collection
        pub_collection = self.generate_collection(self.session, visibility="PUBLIC")
        pub_dataset = self.generate_dataset(
            self.session, collection_visibility="PUBLIC", collection_id=pub_collection.id, published=True
        )
        dep_dir = self.generate_deployment_directory(self.session, pub_dataset.id, upload=True)
        arts = [self.generate_artifact(self.session, pub_dataset.id, ext, upload=False) for ext in DatasetArtifactFileType]

        # Generate revision
        rev_collection = self.generate_collection(self.session, visibility="PRIVATE", id=pub_collection.id)
        rev_datasets = self.generate_dataset(
                self.session,
                collection_visibility="PRIVATE",
                collection_id=rev_collection.id,
                published=True,
                original_id=pub_dataset.id,
            )
        for child in [*arts, dep_dir]:
            self.session.add(clone(child, dataset_id=rev_datasets.id))
        self.session.commit()


        test_url = f"/dp/v1/collections/{pub_collection.id}?visibility=PRIVATE"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        self.app.delete(test_url, headers=headers)
        for dataset in pub_dataset:
            with self.subTest("published artifacts"):
                for artifact in dataset.artifacts:
                    self.assertEqual(self.bucket.Object(pub_dataset.id).get_available_subresources(), 3)
            with self.subTest("published deployed directories"):
                s3_file=pub_dataset.get_cxg_bucket_path(pub_dataset.deploymemnt_directories[0])
                self.assertGreater(self.cellxgene_bucket.Object(s3_file).content_length, 1)
