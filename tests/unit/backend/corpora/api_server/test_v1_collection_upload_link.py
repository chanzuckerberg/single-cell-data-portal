import json
from furl import furl
from mock import patch

from backend.corpora.common.corpora_orm import CollectionVisibility, ProcessingStatus
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.math_utils import GB
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.corpora.api_server.mock_auth import get_auth_token
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup, fixture_file_path
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestCollectionPostUploadLink(BaseAuthAPITest):
    def setUp(self):
        super().setUp()
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        self.dummy_link = "https://www.dropbox.com/s/12345678901234/test.h5ad?dl=0"

    @patch("backend.corpora.common.upload_sfn.start_upload_sfn")
    def test__link__202(self, mocked):
        with EnvironmentSetup({"CORPORA_CONFIG": fixture_file_path("bogo_config.js")}):
            path = "/dp/v1/collections/test_collection_id/upload-links"
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            body = {"url": self.good_link}

            test_url = furl(path=path)
            response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
            self.assertEqual(202, response.status_code)
            self.assertIn("dataset_uuid", json.loads(response.data).keys())

    def test__link_no_auth__401(self):
        path = "/dp/v1/collections/test_collection_id/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json"}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(401, response.status_code)

    @patch(
        "backend.corpora.common.utils.dl_sources.url.DropBoxURL.file_info",
        return_value={"size": 1, "name": "file.h5ad"},
    )
    def test__link_not_owner__403(self, mock_get_file_info):
        path = "/dp/v1/collections/test_collection_id_not_owner/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    def test__bad_link__400(self):
        path = "/dp/v1/collections/test_collection_id/upload-links"

        with self.subTest("Unsupported Provider"):
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            body = {"url": "https://test_url.com"}
            test_url = furl(path=path)
            response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
            self.assertEqual(400, response.status_code)
            self.assertEqual("The dropbox shared link is invalid.", json.loads(response.data)["detail"])

        with self.subTest("Bad Dropbox link"):
            headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
            body = {"url": self.dummy_link}
            test_url = furl(path=path)
            response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
            self.assertEqual(400, response.status_code)
            self.assertEqual("The URL provided causes an error with Dropbox.", json.loads(response.data)["detail"])

    @patch(
        "backend.corpora.common.utils.dl_sources.url.DropBoxURL.file_info", return_value={"size": 1, "name": "file.txt"}
    )
    def test__unsupported_format__400(self, mock_func):
        path = "/dp/v1/collections/test_collection_id/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(400, response.status_code)

    @patch(
        "backend.corpora.common.utils.dl_sources.url.DropBoxURL.file_info",
        return_value={"size": 31 * GB, "name": "file.txt"},
    )
    def test__oversized__413(self, mock_func):
        path = "/dp/v1/collections/test_collection_id/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.dummy_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(413, response.status_code)

    def test__link_fake_collection__403(self):
        path = "/dp/v1/collections/fake/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.good_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)

    def test_link_public_collection__403(self):
        path = "/dp/v1/collections/test_collection_id_public/upload-links"
        headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}
        body = {"url": self.good_link}

        test_url = furl(path=path)
        response = self.app.post(test_url.url, headers=headers, data=json.dumps(body))
        self.assertEqual(403, response.status_code)


class TestCollectionPutUploadLink(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        super().setUp()
        self.good_link = "https://www.dropbox.com/s/ow84zm4h0wkl409/test.h5ad?dl=0"
        self.headers = {"host": "localhost", "Content-Type": "application/json", "Cookie": get_auth_token(self.app)}

    @patch("corpora.common.upload_sfn.start_upload_sfn")
    def test__reupload_published_dataset_during_revision__202(self, mocked):
        """reupload a published dataset during a revision"""
        pub_collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name)
        pub_dataset = self.generate_dataset_with_s3_resources(
            self.session,
            collection_id=pub_collection.id,
            collection_visibility=pub_collection.visibility,
            published=True,
            processing_status={"processing_status": ProcessingStatus.SUCCESS},
        )
        pub_s3_objects = self.get_s3_object_paths_from_dataset(pub_dataset)
        rev_collection = pub_collection.revision()
        path = f"/dp/v1/collections/{pub_collection.id}/upload-links"
        body = {"url": self.good_link, "id": rev_collection.datasets[0].id}

        with EnvironmentSetup({"CORPORA_CONFIG": fixture_file_path("bogo_config.js")}):
            response = self.app.put(path, headers=self.headers, data=json.dumps(body))
            self.assertEqual(202, response.status_code)
            new_datset_id = json.loads(response.body)["dataset_uuid"]
            self.assertEqual(pub_dataset.id, Dataset.get(self.session, new_datset_id).original_id)
            for s3_object in pub_s3_objects:
                self.assertS3FileExists(*s3_object)

    @patch("corpora.common.upload_sfn.start_upload_sfn")
    def test__reupload_unpublished_dataset__202(self, mocked):
        """reupload a unpublished dataset, this removes the old s3 assets. A new uuid is generated"""
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset_with_s3_resources(
            self.session,
            collection_id=collection.id,
            collection_visibility=collection.visibility,
            published=False,
            processing_status={"processing_status": ProcessingStatus.SUCCESS},
        )
        dataset_id = dataset.id
        s3_objects = self.get_s3_object_paths_from_dataset(dataset)
        path = f"/dp/v1/collections/{collection.id}/upload-links"
        body = {"url": self.good_link, "id": dataset_id}

        with EnvironmentSetup({"CORPORA_CONFIG": fixture_file_path("bogo_config.js")}):
            response = self.app.put(path, headers=self.headers, data=json.dumps(body))
            self.assertEqual(202, response.status_code)
            actual_body = json.loads(response.body)
            self.assertEqual(dataset_id, actual_body["dataset_uuid"])
            for s3_object in s3_objects:
                self.assertS3FileDoesNotExist(*s3_object)

    @patch("corpora.common.upload_sfn.start_upload_sfn")
    def test__reupload_unpublished_dataset_during_revision_202(self, mock):
        """reupload a unpublished dataset during a revision, this removes the old s3 assets. A new uuid is generated"""
        collection = self.generate_collection(self.session)
        dataset = self.generate_dataset_with_s3_resources(
            self.session,
            collection_id=collection.id,
            collection_visibility=collection.visibility,
            published=True,
            processing_status={"processing_status": ProcessingStatus.SUCCESS},
        )
        pub_s3_objects = self.get_s3_object_paths_from_dataset(dataset)
        path = f"/dp/v1/collections/{collection.id}/upload-links"
        body = {"url": self.good_link, "id": dataset.id}

        with EnvironmentSetup({"CORPORA_CONFIG": fixture_file_path("bogo_config.js")}):
            response = self.app.put(path, headers=self.headers, data=json.dumps(body))
            self.assertEqual(202, response.status_code)
            for s3_object in pub_s3_objects:
                self.assertS3FileExists(*s3_object)

    def test__reupload_public_dataset__403(self):
        """cannot reupload a public published dataset"""
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PUBLIC.name)
        pub_dataset = self.generate_dataset_with_s3_resources(
            self.session,
            collection_id=collection.id,
            collection_visibility=collection.visibility,
            published=True,
            processing_status={"processing_status": ProcessingStatus.SUCCESS},
        )
        public_dataset_id = pub_dataset.id
        public_s3_objects = self.get_s3_object_paths_from_dataset(pub_dataset)
        path = f"/dp/v1/collections/{collection.id}/upload-links"
        body = {"url": self.good_link, "id": pub_dataset.id}

        with EnvironmentSetup({"CORPORA_CONFIG": fixture_file_path("bogo_config.js")}):
            response = self.app.put(path, headers=self.headers, data=json.dumps(body))
            self.assertEqual(403, response.status_code)
            self.assertIsNotNone(Dataset.get(self.session, public_dataset_id))
            for s3_object in public_s3_objects:
                self.assertS3FileExists(*s3_object)

    def test__reupload_while_processing_dataset__405(self):
        """cannot reupload a dataset that is pending"""
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset_with_s3_resources(
            self.session,
            collection_id=collection.id,
            collection_visibility=collection.visibility,
            processing_status={"processing_status": ProcessingStatus.PENDING},
        )
        s3_objects = self.get_s3_object_paths_from_dataset(dataset)
        path = f"/dp/v1/collections/{collection.id}/upload-links"
        body = {"url": self.good_link, "id": dataset.id}
        with EnvironmentSetup({"CORPORA_CONFIG": fixture_file_path("bogo_config.js")}):
            response = self.app.put(path, headers=self.headers, data=json.dumps(body))
            self.assertEqual(405, response.status_code)
            for s3_object in s3_objects:
                self.assertS3FileExists(*s3_object)

    def test__reupload_dataset_not_owner__403(self):
        collection = self.generate_collection(
            self.session, visibility=CollectionVisibility.PRIVATE.name, owner="someone else"
        )
        dataset = self.generate_dataset_with_s3_resources(
            self.session,
            collection_id=collection.id,
            collection_visibility=collection.visibility,
            published=False,
            processing_status={"processing_status": ProcessingStatus.SUCCESS},
        )
        dataset_id = dataset.id
        path = f"/dp/v1/collections/{collection.id}/upload-links"
        body = {"url": self.good_link, "id": dataset_id}

        with EnvironmentSetup({"CORPORA_CONFIG": fixture_file_path("bogo_config.js")}):
            response = self.app.put(path, headers=self.headers, data=json.dumps(body))
            self.assertEqual(403, response.status_code)

    def test__dataset_not_in_collection__404(self):
        collection = self.generate_collection(self.session, visibility=CollectionVisibility.PRIVATE.name)
        dataset = self.generate_dataset_with_s3_resources(
            self.session,
            published=False,
            processing_status={"processing_status": ProcessingStatus.SUCCESS},
        )
        dataset_id = dataset.id
        path = f"/dp/v1/collections/{collection.id}/upload-links"
        body = {"url": self.good_link, "id": dataset_id}

        with EnvironmentSetup({"CORPORA_CONFIG": fixture_file_path("bogo_config.js")}):
            response = self.app.put(path, headers=self.headers, data=json.dumps(body))
            self.assertEqual(404, response.status_code)
