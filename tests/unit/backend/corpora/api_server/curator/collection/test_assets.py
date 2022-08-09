from unittest.mock import patch, Mock

from backend.corpora.common.corpora_orm import CollectionVisibility, DatasetArtifactFileType
from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.fixtures.config import fake_s3_file
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


# TODO test with dataset_id = test_curator_tag


class TestAsset(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()
        self.test_dataset_id = ["test_dataset_id", "test_curator_tag"]
        self.curator_tag = "curator_tag"
        self.generate_collection(
            self.session,
            id="test_curator_tag_collection_id",
            visibility=CollectionVisibility.PUBLIC.name,
            owner="test_user_id",
            name="test_collection_name",
            description="test_description",
            data_submission_policy_version="0",
            contact_name="Some Body",
            contact_email="somebody@chanzuckerberg.com",
        )
        self.generate_dataset(
            self.session,
            id="test_curator_tag",
            curator_tag="curator_tag",
            revision=0,
            name="test_dataset_name",
            schema_version="2.0.0",
            collection_id="test_curator_tag_collection_id",
            artifacts=[
                dict(
                    filename="test_filename",
                    filetype=DatasetArtifactFileType.H5AD.name,
                    user_submitted=True,
                    s3_uri=fake_s3_file,
                )
            ],
        )

    def test__get_dataset_asset__OK(self):
        bucket = self.CORPORA_TEST_CONFIG["bucket_name"]
        s3_file_name = "test_s3_uri.h5ad"
        content = "Hello world!"
        self.create_s3_object(s3_file_name, bucket, content=content)

        for dataset_id in self.test_dataset_id:
            with self.subTest(dataset_id):
                expected_body = dict(
                    dataset_id=dataset_id,
                    assets=[dict(file_name="test_filename", file_size=len(content), file_type="H5AD")],
                )
                if dataset_id == "test_curator_tag":
                    expected_body["curator_tag"] = self.curator_tag
                    response = self.app.get(
                        "/curation/v1/collections/test_curator_tag_collection_id/datasets/assets",
                        query_string=dict(curator_tag=self.curator_tag),
                    )
                else:
                    response = self.app.get(
                        "/curation/v1/collections/test_collection_id/datasets/assets",
                        query_string=dict(dataset_id=dataset_id),
                    )
                self.assertEqual(200, response.status_code)
                actual_body = response.json
                presign_url = actual_body["assets"][0].pop("presigned_url")
                self.assertIsNotNone(presign_url)
                self.assertEqual(expected_body, actual_body)

    def test__get_dataset_asset__file_error(self):
        for dataset_id in self.test_dataset_id:
            with self.subTest(dataset_id):
                expected_body = dict(
                    dataset_id=dataset_id, assets=[dict(file_name="test_filename", file_size=-1, file_type="H5AD")]
                )
                if dataset_id == "test_curator_tag":
                    expected_body["curator_tag"] = self.curator_tag
                    response = self.app.get(
                        "/curation/v1/collections/test_curator_tag_collection_id/datasets/assets",
                        query_string=dict(curator_tag=self.curator_tag),
                    )
                else:
                    response = self.app.get(
                        "/curation/v1/collections/test_collection_id/datasets/assets",
                        query_string=dict(dataset_id=dataset_id),
                    )
                self.assertEqual(202, response.status_code)
                actual_body = response.json
                presign_url = actual_body["assets"][0].pop("presigned_url", None)
                self.assertIsNone(presign_url)
                self.assertEqual(expected_body, actual_body)

    def test__get_dataset_asset__dataset_NOT_FOUND(self):
        id_or_tag = "bad_id"

        for i in ["uuid", "tag"]:
            with self.subTest(i):
                if i == "tag":
                    response = self.app.get(
                        "/curation/v1/collections/test_curator_tag_collection_id/datasets/assets",
                        query_string=dict(curator_tag=id_or_tag),
                    )
                else:
                    response = self.app.get(
                        "/curation/v1/collections/test_collection_id/datasets/assets",
                        query_string=dict(dataset_id=id_or_tag),
                    )
                self.assertEqual(404, response.status_code)
                actual_body = response.json
                self.assertEqual("Dataset not found.", actual_body["detail"])

    def test__get_dataset_asset__collection_dataset_NOT_FOUND(self):
        """Return Not found when the dataset is not part of the collection"""
        id_or_tag = "bad_id"
        test_url = "/curation/v1/collections/test_collection_id/datasets/assets"
        for i in ["dataset_id", "curator_tag"]:
            with self.subTest(i):
                query = {i: id_or_tag}
                response = self.app.get(test_url, query_string=query)
                self.assertEqual(404, response.status_code)
                actual_body = response.json
                self.assertEqual("Dataset not found.", actual_body["detail"])

    @patch(
        "backend.corpora.lambdas.api.v1.curation.collections.collection_id.assets.get_dataset_else_error",
        return_value=None,
    )
    def test__get_dataset_asset__asset_NOT_FOUND(self, get_dataset_else_error: Mock):
        mocked_dataset = Mock()
        mocked_dataset.get_assets.return_value = None
        get_dataset_else_error.return_value = mocked_dataset
        for dataset_id in self.test_dataset_id:
            with self.subTest(dataset_id):
                if dataset_id == "test_curator_tag":
                    response = self.app.get(
                        "/curation/v1/collections/test_curator_tag_collection_id/datasets/assets",
                        query_string=dict(curator_tag="curator_tag"),
                    )
                else:
                    response = self.app.get(
                        "/curation/v1/collections/test_collection_id/datasets/assets",
                        query_string=dict(dataset_id=dataset_id),
                    )
                self.assertEqual(404, response.status_code)
                actual_body = response.json
                self.assertEqual("No assets found. The dataset may still be processing.", actual_body["detail"])
                mocked_dataset.get_assets.assert_called()
