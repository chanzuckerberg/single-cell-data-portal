from unittest.mock import patch, Mock

from tests.unit.backend.corpora.api_server.base_api_test import BaseAuthAPITest
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS


# TODO test with dataset_uuid = test_curator_tag


class TestAsset(BaseAuthAPITest, CorporaTestCaseUsingMockAWS):
    def setUp(self):
        # Needed for proper setUp resolution in multiple inheritance
        super().setUp()
        self.test_dataset_uuid = ["test_dataset_id", "test_curator_tag"]
        self.curator_tag = "curator_tag"

    def tearDown(self):
        super().tearDown()

    def test__post_dataset_asset__OK(self):
        bucket = self.CORPORA_TEST_CONFIG["bucket_name"]
        s3_file_name = "test_s3_uri.h5ad"
        content = "Hello world!"
        self.create_s3_object(s3_file_name, bucket, content=content)

        for dataset_uuid in self.test_dataset_uuid:
            with self.subTest(dataset_uuid):
                expected_body = dict(
                    dataset_uuid=dataset_uuid,
                    assets=[dict(file_name="test_filename", file_size=len(content), file_type="H5AD")],
                )
                if dataset_uuid == "test_curator_tag":
                    expected_body["curator_tag"] = self.curator_tag
                    response = self.app.post(
                        "/curation/v1/collections/test_collection_id/assets",
                        query_string=dict(curator_tag=self.curator_tag),
                    )
                else:
                    response = self.app.post(
                        "/curation/v1/collections/test_collection_id/assets",
                        query_string=dict(dataset_uuid=dataset_uuid),
                    )
                self.assertEqual(200, response.status_code)
                actual_body = response.json
                presign_url = actual_body["assets"][0].pop("presigned_url")
                self.assertIsNotNone(presign_url)
                self.assertEqual(expected_body, actual_body)

    def test__post_dataset_asset__file_error(self):
        for dataset_uuid in self.test_dataset_uuid:
            with self.subTest(dataset_uuid):
                expected_body = dict(
                    dataset_uuid=dataset_uuid, assets=[dict(file_name="test_filename", file_size=-1, file_type="H5AD")]
                )
                if dataset_uuid == "test_curator_tag":
                    expected_body["curator_tag"] = self.curator_tag
                    response = self.app.post(
                        "/curation/v1/collections/test_collection_id/assets",
                        query_string=dict(curator_tag=self.curator_tag),
                    )
                else:
                    response = self.app.post(
                        "/curation/v1/collections/test_collection_id/assets",
                        query_string=dict(dataset_uuid=dataset_uuid),
                    )
                self.assertEqual(202, response.status_code)
                actual_body = response.json
                presign_url = actual_body["assets"][0].pop("presigned_url", None)
                self.assertIsNone(presign_url)
                self.assertEqual(expected_body, actual_body)

    def test__post_dataset_asset__dataset_NOT_FOUND(self):
        id_or_tag = "bad_id"

        for i in ["uuid", "tag"]:
            with self.subTest(i):
                if i == "tag":
                    response = self.app.post(
                        "/curation/v1/collections/test_collection_id/assets", query_string=dict(curator_tag=id_or_tag)
                    )
                else:
                    response = self.app.post(
                        "/curation/v1/collections/test_collection_id/assets", query_string=dict(dataset_uuid=id_or_tag)
                    )
                self.assertEqual(404, response.status_code)
                actual_body = response.json
                self.assertEqual("Dataset not found.", actual_body["detail"])

    @patch(
        "backend.corpora.lambdas.api.v1.curation.collections.collection_uuid.assets.Dataset.get_assets",
        return_value=None,
    )
    def test__post_dataset_asset__asset_NOT_FOUND(self, get_asset: Mock):
        for dataset_uuid in self.test_dataset_uuid:
            with self.subTest(dataset_uuid):
                if dataset_uuid == "test_curator_tag":
                    response = self.app.post(
                        "/curation/v1/collections/test_collection_id/assets",
                        query_string=dict(curator_tag="curator_tag"),
                    )
                else:
                    response = self.app.post(
                        "/curation/v1/collections/test_collection_id/assets",
                        query_string=dict(dataset_uuid=dataset_uuid),
                    )
                self.assertEqual(404, response.status_code)
                actual_body = response.json
                self.assertEqual(
                    "No assets found not found. The dataset may still be processing.", actual_body["detail"]
                )
                get_asset.assert_called()
