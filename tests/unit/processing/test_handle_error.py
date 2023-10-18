import os
from typing import Dict
from unittest.mock import Mock, patch

import pytest

from backend.layers.common.entities import CollectionId, CollectionVersionId, DatasetStatus, DatasetVersionId
from backend.layers.processing.upload_failures.app import (
    FAILED_ARTIFACT_CLEANUP_MESSAGE,
    FAILED_CXG_CLEANUP_MESSAGE,
    FAILED_DATASET_CLEANUP_MESSAGE,
    cleanup_artifacts,
    get_failure_slack_notification_message,
    parse_event,
)

module_path = "backend.layers.processing.upload_failures.app"


@pytest.fixture
def sample_slack_header_block():
    return {
        "type": "header",
        "text": {
            "type": "plain_text",
            "text": "Dataset failed to process:fire:",
            "emoji": True,
        },
    }


@pytest.fixture
def sample_slack_status_block():
    return {
        "type": "section",
        "text": {
            "type": "mrkdwn",
            "text": "```{\n"
            '  "cxg_status": null,\n'
            '  "h5ad_status": null,\n'
            '  "processing_status": null,\n'
            '  "rds_status": null,\n'
            '  "upload_status": null,\n'
            '  "validation_message": null,\n'
            '  "validation_status": null\n'
            "}```",
        },
    }


@pytest.fixture
def sample_slack_status_block_empty():
    return {
        "type": "section",
        "text": {
            "type": "mrkdwn",
            "text": "``````",
        },
    }


def test_parse_event_with_empty_event():
    (
        dataset_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event({})

    assert dataset_id is None
    assert collection_version_id is None
    assert error_step_name is None
    assert error_job_id is None
    assert error_aws_regions is None
    assert error_cause == ""
    assert execution_arn is None


def test_parse_event_with_error_cause():
    expected_error_cause = (
        '{"JobName": "Step1", "JobId": "789", "Container": {"Environment": [{"Name": "AWS_DEFAULT_REGION", '
        '"Value": "us-east-1"}]}}'
    )
    event = {
        "execution_id": "arn",
        "dataset_id": "123",
        "collection_id": "456",
        "error": {"Cause": expected_error_cause},
    }
    (
        dataset_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event(event)

    assert dataset_id == "123"
    assert collection_version_id == "456"
    assert error_step_name == "Step1"
    assert error_job_id == "789"
    assert error_aws_regions == "us-east-1"
    assert expected_error_cause is error_cause
    assert execution_arn == "arn"


def test_parse_event_without_error_cause():
    event = {"dataset_id": "123", "collection_id": "456", "error": {}}

    (
        dataset_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event(event)

    assert dataset_id == "123"
    assert collection_version_id == "456"
    assert error_step_name is None
    assert error_job_id is None
    assert error_aws_regions is None
    assert error_cause == ""
    assert execution_arn is None


def test_parse_event_with_invalid_error_cause():
    event = {"dataset_id": "123", "collection_id": "456", "error": {"Cause": "invalid JSON"}}

    (
        dataset_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event(event)

    assert dataset_id == "123"
    assert collection_version_id == "456"
    assert error_step_name is None
    assert error_job_id is None
    assert error_aws_regions is None
    assert error_cause == "invalid JSON"
    assert execution_arn is None


def mock_get_dataset_version(collection_id):
    MockDatasetVersionId = Mock()
    MockDatasetVersionId.collection_id = collection_id
    MockDatasetVersionId.status = DatasetStatus(None, None, None, None, None, None)
    return MockDatasetVersionId


def test_get_failure_slack_notification_message_with_dataset_id_none(
    sample_slack_header_block, sample_slack_status_block_empty
):
    dataset_id = None
    step_name = "Step 1"
    job_id = "123456"
    aws_regions = "us-west-2"
    execution_arn = "arn:aws:states:us-west-2:123456789012:execution:MyStateMachine"
    collection_version_id = "collection_version_id123"

    with patch(f"{module_path}.logger") as logger_mock:
        result = get_failure_slack_notification_message(
            dataset_id, collection_version_id, step_name, job_id, aws_regions, execution_arn
        )
    assert result == {
        "blocks": [
            sample_slack_header_block,
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! @sc-oncall-eng please follow the [triage steps](https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit#bookmark=id.3ofm47y0709y)\n"
                    "*Owner*: \n"
                    f"*Collection Version URL*: https://cellxgene.cziscience.com/collections/{collection_version_id}\n"
                    "*Batch Job ID*: <https://us-west-2.console.aws.amazon.com/batch/v2/home?region=us-west-2"
                    "#jobs/detail/123456|123456>\n"
                    "*Step Function ARN*: "
                    "<https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/v2/executions"
                    "/details/arn:aws:states:us-west-2:123456789012:execution:MyStateMachine|arn:aws:states"
                    ":us-west-2:123456789012:execution:MyStateMachine>\n"
                    f"*Error Step*: {step_name}\n"
                    f"*Dataset ID*: None(not found)\n"
                    f"*Processing Status*:\n",
                },
            },
            sample_slack_status_block_empty,
        ]
    }
    logger_mock.error.assert_called_with("Dataset None not found")


def test_get_failure_slack_notification_message_with_dataset_not_found(
    sample_slack_header_block, sample_slack_status_block_empty
):
    dataset_id = "dataset123"
    step_name = "Step 1"
    job_id = "123456"
    aws_regions = "us-west-2"
    execution_arn = "arn:aws:states:us-west-2:123456789012:execution:MyStateMachine"
    collection_version_id = "collection_version_id123"

    get_dataset_version_mock = Mock(return_value=None)
    get_business_logic_mock = Mock()
    get_business_logic_mock.get_dataset_version = get_dataset_version_mock
    get_business_logic_constructor_mock = Mock(return_value=get_business_logic_mock)

    logger_mock = Mock()

    with patch(f"{module_path}.get_business_logic", get_business_logic_constructor_mock), patch(
        f"{module_path}.logger", logger_mock
    ):
        result = get_failure_slack_notification_message(
            dataset_id, collection_version_id, step_name, job_id, aws_regions, execution_arn
        )

    assert result == {
        "blocks": [
            sample_slack_header_block,
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! @sc-oncall-eng please follow the [triage steps](https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit#bookmark=id.3ofm47y0709y)\n"
                    "*Owner*: \n"
                    f"*Collection Version URL*: https://cellxgene.cziscience.com/collections/{collection_version_id}\n"
                    "*Batch Job ID*: <https://us-west-2.console.aws.amazon.com/batch/v2/home?region=us-west-2"
                    "#jobs/detail/123456|123456>\n"
                    "*Step Function ARN*: "
                    "<https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/v2/executions"
                    "/details/arn:aws:states:us-west-2:123456789012:execution:MyStateMachine|arn:aws:states"
                    ":us-west-2:123456789012:execution:MyStateMachine>\n"
                    f"*Error Step*: {step_name}\n"
                    f"*Dataset ID*: {dataset_id}(not found)\n"
                    f"*Processing Status*:\n",
                },
            },
            sample_slack_status_block_empty,
        ]
    }
    logger_mock.error.assert_called_with("Dataset dataset123 not found")
    get_dataset_version_mock.assert_called_with(DatasetVersionId(dataset_id))


def mock_collection_version(owner, version_id):
    MockCollectionVersion = Mock()
    MockCollectionVersion.owner = owner
    MockCollectionVersion.version_id = version_id
    return MockCollectionVersion


def test_get_failure_slack_notification_message_with_missing_collection(
    sample_slack_header_block, sample_slack_status_block
):
    dataset_id = "dataset123"
    collection_id = "collection123"
    collection_version_id = "collection_version_id123"
    step_name = "Step 1"
    job_id = "123456"
    aws_regions = "us-west-2"
    execution_arn = "arn:aws:states:us-west-2:123456789012:execution:MyStateMachine"

    get_dataset_version_mock = Mock(return_value=mock_get_dataset_version(CollectionId(collection_id)))
    get_unpublished_collection_version_from_canonical_mock = Mock(return_value=None)

    get_business_logic_mock = Mock()
    get_business_logic_mock.get_dataset_version = get_dataset_version_mock
    get_business_logic_mock.get_unpublished_collection_version_from_canonical = (
        get_unpublished_collection_version_from_canonical_mock
    )
    get_business_logic_constructor_mock = Mock(return_value=get_business_logic_mock)

    logger_mock = Mock()

    with patch(f"{module_path}.get_business_logic", get_business_logic_constructor_mock), patch(
        f"{module_path}.logger", logger_mock
    ):
        result = get_failure_slack_notification_message(
            dataset_id, collection_version_id, step_name, job_id, aws_regions, execution_arn
        )

    assert result == {
        "blocks": [
            sample_slack_header_block,
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! @sc-oncall-eng please follow the [triage steps](https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit#bookmark=id.3ofm47y0709y)\n"
                    f"*Owner*: \n"
                    f"*Collection Version URL*: https://cellxgene.cziscience.com/collections/{collection_version_id}\n"
                    "*Batch Job ID*: <https://us-west-2.console.aws.amazon.com/batch/v2/home?region=us-west-2"
                    "#jobs/detail/123456|123456>\n"
                    "*Step Function ARN*: "
                    "<https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/v2/executions"
                    "/details/arn:aws:states:us-west-2:123456789012:execution:MyStateMachine|arn:aws:states"
                    ":us-west-2:123456789012:execution:MyStateMachine>\n"
                    f"*Error Step*: {step_name}\n"
                    f"*Dataset ID*: {dataset_id}\n"
                    f"*Processing Status*:\n",
                },
            },
            sample_slack_status_block,
        ]
    }
    logger_mock.error.assert_called_with(f"Collection {collection_id} not found")
    get_dataset_version_mock.assert_called_with(DatasetVersionId(dataset_id))
    get_unpublished_collection_version_from_canonical_mock.assert_called_with(CollectionId(collection_id))


def test_get_failure_slack_notification_message_with_dataset_and_collection(
    sample_slack_header_block, sample_slack_status_block
):
    dataset_id = "dataset123"
    collection_id = "collection123"
    step_name = "Step 1"
    job_id = "123456"
    aws_regions = "us-west-2"
    execution_arn = "arn:aws:states:us-west-2:123456789012:execution:MyStateMachine"
    owner = "test"
    collection_version_id = "version123"

    get_dataset_version_mock = Mock(return_value=mock_get_dataset_version(CollectionId(collection_id)))
    get_unpublished_collection_version_from_canonical_mock = Mock(
        return_value=mock_collection_version(owner, CollectionVersionId(collection_version_id))
    )

    get_business_logic_mock = Mock()
    get_business_logic_mock.get_dataset_version = get_dataset_version_mock
    get_business_logic_mock.get_unpublished_collection_version_from_canonical = (
        get_unpublished_collection_version_from_canonical_mock
    )
    get_business_logic_constructor_mock = Mock(return_value=get_business_logic_mock)

    with patch(f"{module_path}.get_business_logic", get_business_logic_constructor_mock):
        result = get_failure_slack_notification_message(
            dataset_id, collection_version_id, step_name, job_id, aws_regions, execution_arn
        )

    assert result == {
        "blocks": [
            sample_slack_header_block,
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! @sc-oncall-eng please follow the [triage steps](https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit#bookmark=id.3ofm47y0709y)\n"
                    f"*Owner*: {owner}\n"
                    f"*Collection Version URL*: https://cellxgene.cziscience.com/collections/{collection_version_id}\n"
                    "*Batch Job ID*: <https://us-west-2.console.aws.amazon.com/batch/v2/home?region=us-west-2"
                    "#jobs/detail/123456|123456>\n"
                    "*Step Function ARN*: "
                    "<https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/v2/executions"
                    "/details/arn:aws:states:us-west-2:123456789012:execution:MyStateMachine|arn:aws:states"
                    ":us-west-2:123456789012:execution:MyStateMachine>\n"
                    f"*Error Step*: {step_name}\n"
                    f"*Dataset ID*: {dataset_id}\n"
                    f"*Processing Status*:\n",
                },
            },
            sample_slack_status_block,
        ]
    }
    get_dataset_version_mock.assert_called_with(DatasetVersionId(dataset_id))
    get_unpublished_collection_version_from_canonical_mock.assert_called_with(CollectionId(collection_id))


@pytest.fixture
def mock_env_vars() -> Dict[str, str]:
    mock_env_vars = {
        "ARTIFACT_BUCKET": "artifact_bucket",
        "DATASETS_BUCKET": "datasets_bucket",
        "CELLXGENE_BUCKET": "cxg_bucket",
    }
    with patch.dict(os.environ, mock_env_vars):
        yield mock_env_vars


@pytest.fixture
def mock_delete_many_from_s3() -> Mock:
    with patch(f"{module_path}.delete_many_from_s3") as mock_delete_many_from_s3:
        yield mock_delete_many_from_s3


@pytest.fixture
def dataset_id() -> str:
    return "example_dataset"


class TestCleanupArtifacts:
    @pytest.mark.parametrize("error_step", ["download-validate", "", None])
    def test_cleanup_artifacts__OK(self, mock_env_vars, mock_delete_many_from_s3, dataset_id, error_step):
        """Check that all artifacts are deleted for the given cases."""
        cleanup_artifacts(dataset_id, error_step)

        # Assertions
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["ARTIFACT_BUCKET"], dataset_id + "/")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["DATASETS_BUCKET"], dataset_id + ".")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["CELLXGENE_BUCKET"], dataset_id + ".cxg/")
        assert mock_delete_many_from_s3.call_count == 3

    def test_cleanup_artifacts__not_download_validate(self, mock_env_vars, mock_delete_many_from_s3, dataset_id):
        """Check that file in the artifact bucket are not delete if error_step is not download-validate."""
        cleanup_artifacts(dataset_id, "not_download_validate")

        # Assertions
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["DATASETS_BUCKET"], dataset_id + ".")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["CELLXGENE_BUCKET"], dataset_id + ".cxg/")
        assert mock_delete_many_from_s3.call_count == 2

    @patch.dict(os.environ, clear=True)
    def test_cleanup_artifacts__no_buckets(self, caplog, mock_delete_many_from_s3, dataset_id):
        """Check that no files are deleted if buckets are not specified."""
        cleanup_artifacts(dataset_id)

        # Assertions
        mock_delete_many_from_s3.assert_not_called()
        assert FAILED_ARTIFACT_CLEANUP_MESSAGE in caplog.text
        assert FAILED_CXG_CLEANUP_MESSAGE in caplog.text
        assert FAILED_DATASET_CLEANUP_MESSAGE in caplog.text

    def test_cleanup_artifacts__elete_many_from_s3_error(
        self, caplog, mock_env_vars, mock_delete_many_from_s3, dataset_id
    ):
        """Check that delete_many_from_s3 errors are logged but do not raise exceptions."""
        mock_delete_many_from_s3.side_effect = Exception("Boom!")
        cleanup_artifacts(dataset_id)

        # Assertions
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["ARTIFACT_BUCKET"], dataset_id + "/")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["DATASETS_BUCKET"], dataset_id + ".")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["CELLXGENE_BUCKET"], dataset_id + ".cxg/")
        assert mock_delete_many_from_s3.call_count == 3
        assert FAILED_ARTIFACT_CLEANUP_MESSAGE in caplog.text
        assert FAILED_CXG_CLEANUP_MESSAGE in caplog.text
        assert FAILED_DATASET_CLEANUP_MESSAGE in caplog.text
