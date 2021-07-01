import mock

from botocore.exceptions import ClientError

from backend.corpora.common.corpora_orm import (
    DatasetArtifactType,
    DatasetArtifactFileType,
)
from backend.corpora.common.entities.dataset_asset import DatasetAsset
from backend.corpora.common.utils.db_session import DBSessionMaker
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from tests.unit.backend.fixtures.generate_data_mixin import GenerateDataMixin


class TestDatasetAsset(CorporaTestCaseUsingMockAWS, GenerateDataMixin):
    def setUp(self):
        super().setUp()
        self.uuid = "test_dataset_artifact_id"
        self.bucket_name = self.CORPORA_TEST_CONFIG["bucket_name"]
        self.session = DBSessionMaker().session()

    def tearDown(self):
        self.session.close()
        super().tearDown()

    def test__get__ok(self):
        asset = DatasetAsset.get(self.session, self.uuid)
        self.assertEqual(self.uuid, asset.id)
        self.assertEqual(self.CORPORA_TEST_CONFIG["bucket_name"], asset.bucket_name)
        self.assertEqual("test_s3_uri.h5ad", asset.key_name)

    def test__generate_file_url__OK(self):
        file_name = "test_generate_file_url.h5ad"

        # Create the Dataset Asset
        asset = self.create_dataset_asset(file_name)

        url = asset.generate_file_url()
        self.assertIn(file_name, url)
        self.assertIn("Expires=", url)

    @mock.patch("backend.corpora.common.entities.dataset_asset.s3_client")
    def test__generate_file_url__ERROR(self, mock_s3):
        mock_s3.generate_presigned_url.side_effect = ClientError({}, "mock ClientError")
        # Create the Dataset Asset
        asset = self.create_dataset_asset("test__generate_file_url__ERROR.h5ad")

        with self.assertLogs(level="ERROR") as logger:
            url = asset.generate_file_url()
        self.assertIn("Failed to generate presigned URL", logger.output[0])
        self.assertEqual(None, url)

    def test__get_file_size__OK(self):
        file_name = "test_head_file.txt"
        content = "This is test_head_file."

        # Create S3 object
        self.create_s3_object(file_name, self.bucket_name, content=content)

        # Create the Dataset Asset
        asset = self.create_dataset_asset(file_name)

        file_size = asset.get_file_size()
        self.assertEqual(len(content), file_size)

    def test__get_file_size__not_found(self):
        file_name = "fake_file.txt"

        # Create the Dataset Asset
        asset = self.create_dataset_asset(file_name)

        with self.assertLogs(level="ERROR") as logger:
            file_size = asset.get_file_size()
        self.assertIn("Failed to retrieve meta data", logger.output[0])
        self.assertIsNone(file_size)

    def test__delete_dataset_asset__ok(self):
        file_name = "test_head_file.txt"
        content = "This is test_head_file."

        # Create S3 object
        self.create_s3_object(file_name, self.bucket_name, content=content)

        # Create the Dataset Asset
        asset = self.create_dataset_asset(file_name)

        # Check file exists in s3 by getting the file size
        file_size = asset.get_file_size()
        self.assertEqual(len(content), file_size)

        # Delete the asset
        asset.delete_from_s3()

        # Check file is no longer in s3
        with self.assertLogs(level="ERROR") as logger:
            file_size = asset.get_file_size()
        self.assertIn("Failed to retrieve meta data", logger.output[0])

    def test__delete_dataset_asset__not_found(self):
        file_name = "test_head_file.txt"
        content = "This is test_head_file."

        # Create S3 object
        self.create_s3_object(file_name, self.bucket_name, content=content)

        # Create the Dataset Asset
        asset = self.create_dataset_asset(file_name)

        # Check file exists in s3 by getting the file size
        file_size = asset.get_file_size()
        self.assertEqual(len(content), file_size)

        # Delete the asset
        asset.delete_from_s3()

        # Check file is no longer in s3
        with self.assertLogs(level="ERROR") as logger:
            file_size = asset.get_file_size()
        self.assertIn("Failed to retrieve meta data", logger.output[0])

        # Delete an asset that no longer exists - this should not raise an error
        asset.delete_from_s3()

    def create_dataset_asset(self, file_name):
        # Create the Dataset Asset
        artifact_params = dict(
            filename="filename_1",
            filetype=DatasetArtifactFileType.H5AD,
            type=DatasetArtifactType.ORIGINAL,
            user_submitted=True,
            s3_uri=f"s3://{self.bucket_name}/{file_name}",
        )
        dataset = self.generate_dataset(self.session, artifacts=[artifact_params])
        return dataset.get_asset(dataset.artifacts[0].id)

    def test_get_most_recent_artifact(self):
        first_uri = "some_uri_0"
        second_uri = "some_uri_1"
        artifact_params_0 = dict(
            filename="filename_x",
            filetype=DatasetArtifactFileType.CXG,
            type=DatasetArtifactType.ORIGINAL,
            user_submitted=True,
            s3_uri=first_uri,
        )
        dataset = self.generate_dataset(self.session, artifacts=[artifact_params_0])
        DatasetAsset.create(
            self.session,
            dataset_id=dataset.id,
            filename="some file",
            filetype=DatasetArtifactFileType.CXG,
            type_enum=DatasetArtifactType.REMIX,
            user_submitted=True,
            s3_uri=second_uri,
        )
        artifact = dataset.get_most_recent_artifact(DatasetArtifactFileType.CXG)
        self.assertEqual(artifact.s3_uri, second_uri)
