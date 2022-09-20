from unittest.mock import patch, Mock

from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


class TestAsset(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()
        self.test_dataset_id = "test_dataset_id"

    def test__get_dataset_asset__OK(self):
        bucket = self.CORPORA_TEST_CONFIG["bucket_name"]
        s3_file_name = "test_s3_uri.h5ad"
        content = "Hello world!"
        self.create_s3_object(s3_file_name, bucket, content=content)

        expected_body = dict(
            dataset_id=self.test_dataset_id,
            assets=[dict(filename="test_filename", filesize=len(content), filetype="H5AD")],
        )

        response = self.app.get(
            f"/curation/v1/collections/test_collection_id/datasets/{self.test_dataset_id}/assets",
            query_string=dict(dataset_id=self.test_dataset_id),
        )
        self.assertEqual(200, response.status_code)
        actual_body = response.json
        presign_url = actual_body["assets"][0].pop("presigned_url")
        self.assertIsNotNone(presign_url)
        self.assertEqual(expected_body, actual_body)

    def test__get_dataset_asset__file_error(self):
        expected_body = dict(
            dataset_id=self.test_dataset_id, assets=[dict(filename="test_filename", filesize=-1, filetype="H5AD")]
        )

        response = self.app.get(
            f"/curation/v1/collections/test_collection_id/datasets/{self.test_dataset_id}/assets",
        )
        self.assertEqual(202, response.status_code)
        actual_body = response.json
        presign_url = actual_body["assets"][0].pop("presigned_url", None)
        self.assertIsNone(presign_url)
        self.assertEqual(expected_body, actual_body)

    def test__get_dataset_asset__dataset_NOT_FOUND(self):
        bad_id = "bad_id"

        response = self.app.get(
            f"/curation/v1/collections/test_collection_id/datasets/{bad_id}/assets",
        )
        self.assertEqual(404, response.status_code)
        actual_body = response.json
        self.assertEqual("Dataset not found.", actual_body["detail"])

    def test__get_dataset_asset__collection_dataset_NOT_FOUND(self):
        """Return Not found when the dataset is not part of the collection"""
        bad_id = "bad_id"
        test_url = f"/curation/v1/collections/test_collection_id/datasets/{bad_id}/assets"
        response = self.app.get(test_url)
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
        response = self.app.get(f"/curation/v1/collections/test_collection_id/datasets/{self.test_dataset_id}/assets")
        self.assertEqual(404, response.status_code)
        actual_body = response.json
        self.assertEqual("No assets found. The dataset may still be processing.", actual_body["detail"])
        mocked_dataset.get_assets.assert_called()
