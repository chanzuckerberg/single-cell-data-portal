import json
from unittest.mock import MagicMock, Mock, patch

import pytest
import scanpy

from backend.common.utils.math_utils import GB
from backend.layers.common.entities import DatasetArtifactType, DatasetUploadStatus
from backend.layers.processing.process_download import ProcessDownload
from tests.unit.processing.base_processing_test import BaseProcessingTest

test_environment = {"REMOTE_DEV_PREFIX": "fake-stack", "DEPLOYMENT_STAGE": "test"}


class TestProcessDownload(BaseProcessingTest):
    @patch("os.environ", test_environment)
    @patch("backend.layers.processing.process_download.StepFunctionProvider")
    @patch("scanpy.read_h5ad")
    def test_process_download_success(self, mock_read_h5ad, mock_sfn_provider):
        """
        ProcessValidate should:
        1. Download the h5ad artifact
        2. Set upload status to UPLOADED
        3. upload the original file to S3
        """
        dropbox_uri = "https://www.dropbox.com/s/fake_location/test.h5ad?dl=0"
        bucket_name = "fake_bucket_name"
        stack_name = test_environment["REMOTE_DEV_PREFIX"]
        deployment_stage = test_environment["DEPLOYMENT_STAGE"]
        collection = self.generate_unpublished_collection()
        dataset_vid, dataset_id = self.business_logic.ingest_dataset(collection.version_id, dropbox_uri, None, None)
        # Mock anndata object
        mock_anndata = Mock(spec=scanpy.AnnData)
        mock_anndata.n_obs = 10000
        mock_anndata.n_vars = 10000
        mock_read_h5ad.return_value = mock_anndata

        # Mock SFN client
        mock_sfn = Mock()
        mock_sfn_provider.return_value = mock_sfn

        # This is where we're at when we start the SFN
        pdv = ProcessDownload(self.business_logic, self.uri_provider, self.s3_provider)
        pdv.process(dataset_vid, dropbox_uri, bucket_name, "fake_sfn_task_token")

        status = self.business_logic.get_dataset_status(dataset_vid)
        self.assertEqual(status.upload_status, DatasetUploadStatus.UPLOADED)

        # Assert mocks
        mock_read_h5ad.assert_called_with("raw.h5ad")
        mock_sfn.client.send_task_success.assert_called_with(
            taskToken="fake_sfn_task_token",
            output=json.dumps(
                {
                    "JobDefinitionName": f"dp-{deployment_stage}-{stack_name}-ingest-process-{dataset_vid.id}",
                    "Vcpus": 1,
                    "Memory": 4000,
                    "LinuxParameters": {"Swappiness": 60, "MaxSwap": 0},
                }
            ),
        )

        # Verify that both the original (raw.h5ad) and the labeled (local.h5ad) files are there
        self.assertTrue(self.s3_provider.uri_exists(f"s3://{bucket_name}/{stack_name}/{dataset_vid.id}/raw.h5ad"))

        artifacts = list(self.business_logic.get_dataset_artifacts(dataset_vid))
        self.assertEqual(1, len(artifacts))
        artifact = artifacts[0]
        artifact.type = DatasetArtifactType.RAW_H5AD
        artifact.uri = f"s3://fake_bucket_name/{stack_name}/{dataset_vid.id}/raw.h5ad"

    def test_download_from_s3_uri(self):
        """
        Call process download using an s3 uri
        """

        s3_uri = "s3://fake_bucket_name/fake_key/fake_file.h5ad"
        pdv = ProcessDownload(Mock(), self.uri_provider, Mock())
        pdv.download_from_s3 = Mock()

        assert pdv.download_from_source_uri(s3_uri, "fake_local_path") == "fake_local_path"
        pdv.download_from_s3.assert_called_once_with(
            bucket_name="fake_bucket_name", object_key="fake_key/fake_file.h5ad", local_filename="fake_local_path"
        )

    def test_download_from_dropbox_uri(self):
        """
        Call process download using a dropbox uri
        """

        dropbox_uri = "https://www.dropbox.com/s/fake_location/test.h5ad?dl=1"
        pdv = ProcessDownload(Mock(), self.uri_provider, Mock())
        pdv.download = Mock()

        assert pdv.download_from_source_uri(dropbox_uri, "fake_local_path") == "fake_local_path"
        pdv.download.assert_called_once_with(dropbox_uri, "fake_local_path")

    def test_download_unknown_uri(self):
        """
        Call process download using unknown
        """

        uri = "fake://fake_bucket_name/fake_key/fake_file.h5ad"
        pdv = ProcessDownload(Mock(), self.uri_provider, Mock())
        pdv.download_from_s3 = Mock()
        with pytest.raises(ValueError, match=f"Malformed source URI: {uri}"):
            pdv.download_from_source_uri(uri, "fake_local_path")


@pytest.fixture
def mock_ProcessDownload():
    return ProcessDownload(Mock(), Mock(), Mock())


def sample_adata(n_obs: int, n_vars: int):
    # Create a sample AnnData object for testing
    adata = MagicMock(spec=scanpy.AnnData, n_obs=n_obs, n_vars=n_vars)
    return adata


@pytest.fixture
def mock_read_h5ad():
    with patch("scanpy.read_h5ad") as mock_read_h5ad:
        mock_read_h5ad.return_value = sample_adata(1, 2 * GB)
        yield mock_read_h5ad


# Arrange
@pytest.mark.parametrize(
    "adata, memory_modifier, expected",
    [
        (sample_adata(1, 2 * GB), 1, {"Vcpus": 1, "Memory": 4000, "MaxSwap": 0}),  # minimum memory
        (sample_adata(1, 5 * GB), 1, {"Vcpus": 2, "Memory": 5120, "MaxSwap": 0}),  # above minimum memory
        (sample_adata(1, 5 * GB), 1.5, {"Vcpus": 2, "Memory": 7680, "MaxSwap": 0}),  # modifier adjusted
        (sample_adata(1, 64 * GB), 1, {"Vcpus": 16, "Memory": 64000, "MaxSwap": 300000}),  # maximum memory
    ],
)
def test_estimate_resource_requirements_positive(mock_ProcessDownload, adata, memory_modifier, expected):
    # Act & Assert
    assert expected == mock_ProcessDownload.estimate_resource_requirements(adata, memory_modifier=memory_modifier)


@pytest.mark.parametrize(
    "environ,expected",
    [
        ({"REMOTE_DEV_PREFIX": "/stack/", "DEPLOYMENT_STAGE": "test"}, "dp-test-stack-ingest-process-fake_dataset_id"),
        ({"REMOTE_DEV_PREFIX": "stack", "DEPLOYMENT_STAGE": "test"}, "dp-test-stack-ingest-process-fake_dataset_id"),
        ({"DEPLOYMENT_STAGE": "test"}, "dp-test-ingest-process-fake_dataset_id"),
    ],
)
def test_get_job_definion_name(mock_ProcessDownload, environ, expected):
    # Arrange
    with patch("os.environ", environ):
        dataset_id = "fake_dataset_id"

        # Act
        result = mock_ProcessDownload.get_job_definion_name(dataset_id)

        # Assert
        assert result == expected


def test_remove_prefix(mock_ProcessDownload):
    # Act & Assert
    assert mock_ProcessDownload.remove_prefix("prefixfake", "prefix") == "fake"


def test_create_batch_job_definition_parameters(mock_ProcessDownload, mock_read_h5ad):
    # Arrange
    mock_ProcessDownload.get_job_definion_name = Mock(return_value="fake_job_definition_name")
    mock_ProcessDownload.estimate_resource_requirements = Mock(return_value={"Vcpus": 1, "Memory": 4000, "MaxSwap": 0})

    # Act
    resp = mock_ProcessDownload.create_batch_job_definition_parameters("local_file.h5ad", "fake_dataset_id")

    # Assert
    mock_ProcessDownload.estimate_resource_requirements.assert_called_once_with(mock_read_h5ad.return_value)
    mock_ProcessDownload.get_job_definion_name.assert_called_once_with("fake_dataset_id")
    assert resp == {
        "JobDefinitionName": "fake_job_definition_name",
        "Vcpus": 1,
        "Memory": 4000,
        "LinuxParameters": {
            "Swappiness": 60,
            "MaxSwap": 0,
        },
    }
