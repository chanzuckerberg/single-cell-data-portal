import datetime
import logging
import os
from typing import Dict
from unittest.mock import Mock, patch

import pytest

from backend.layers.common.entities import (
    CanonicalCollection,
    CollectionId,
    CollectionMetadata,
    CollectionVersionId,
    CollectionVersionWithDatasets,
    DatasetArtifact,
    DatasetArtifactId,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetStatus,
    DatasetVersionId,
)
from backend.layers.processing.upload_failures.app import (
    FAILED_ARTIFACT_CLEANUP_MESSAGE,
    FAILED_ATAC_DATASET_MESSAGE,
    FAILED_CXG_CLEANUP_MESSAGE,
    FAILED_DATASET_CLEANUP_MESSAGE,
    cleanup_artifacts,
    delete_atac_fragment_files,
    get_failure_slack_notification_message,
    handle_failure,
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
            '  "atac_status": null,\n'
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


@pytest.fixture
def get_collection_version_mock():
    return Mock(
        return_value=CollectionVersionWithDatasets(
            datasets=[],
            collection_id=CollectionId("collection123"),
            version_id=CollectionVersionId("collection_version_id123"),
            owner="mock owner",
            curator_name="mock curator",
            metadata=CollectionMetadata(
                name="mock name",
                description="mock description",
                contact_name="mock contact name",
                contact_email="mock contact email",
                links=[],
            ),
            publisher_metadata={},
            published_at=None,
            created_at=datetime.datetime.now(),
            schema_version="5.1.0",
            canonical_collection=CanonicalCollection(
                id=CollectionId("collection123"),
                version_id=None,
                originally_published_at=None,
                revised_at=None,
                tombstoned=False,
            ),
            has_custom_dataset_order=False,
            is_auto_version=False,
            data_submission_policy_version="2.0",
        )
    )


@pytest.fixture
def mock_business_logic(monkeypatch):
    business_logic_mock = Mock()
    monkeypatch.setattr(f"{module_path}.get_business_logic", Mock(return_value=business_logic_mock))
    return business_logic_mock


def test_parse_event_with_empty_event():
    (
        dataset_version_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event({})

    assert dataset_version_id is None
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
        "dataset_version_id": "123",
        "collection_id": "789",
        "collection_version_id": "456",
        "error": {"Cause": expected_error_cause},
    }
    (
        dataset_version_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event(event)

    assert dataset_version_id == "123"
    assert collection_version_id == "456"
    assert error_step_name == "Step1"
    assert error_job_id == "789"
    assert error_aws_regions == "us-east-1"
    assert expected_error_cause is error_cause
    assert execution_arn == "arn"


def test_parse_event_without_error_cause():
    event = {"dataset_version_id": "123", "collection_version_id": "456", "error": {}}

    (
        dataset_version_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event(event)

    assert dataset_version_id == "123"
    assert collection_version_id == "456"
    assert error_step_name is None
    assert error_job_id is None
    assert error_aws_regions is None
    assert error_cause == ""
    assert execution_arn is None


def test_parse_event_with_invalid_error_cause():
    event = {
        "dataset_version_id": "123",
        "collection_version_id": "456",
        "error": {"Cause": "invalid JSON"},
    }

    (
        dataset_version_id,
        collection_version_id,
        error_step_name,
        error_job_id,
        error_aws_regions,
        error_cause,
        execution_arn,
    ) = parse_event(event)

    assert dataset_version_id == "123"
    assert collection_version_id == "456"
    assert error_step_name is None
    assert error_job_id is None
    assert error_aws_regions is None
    assert error_cause == "invalid JSON"
    assert execution_arn is None


def mock_get_dataset_version(collection_id):
    MockDatasetVersionId = Mock()
    MockDatasetVersionId.collection_id = collection_id
    MockDatasetVersionId.status = DatasetStatus(*[None] * 7)
    return MockDatasetVersionId


def test_migration_event_does_not_trigger_slack(mock_business_logic):
    mock_trigger_slack = Mock()
    mock_context = Mock()
    with patch("backend.layers.processing.upload_failures.app.trigger_slack_notification", mock_trigger_slack):
        event = {
            "dataset_version_id": "123",
            "collection_version_id": "456",
            "error": {},
            "execution_id": "arn:aws:states:us-west-2:migrate_123456789012:execution:MyStateMachine",
        }
        handle_failure(event, mock_context, delete_artifacts=False)
        mock_trigger_slack.assert_not_called()


def test_non_migration_event_triggers_slack(mock_business_logic):
    mock_trigger_slack = Mock()
    mock_context = Mock()
    with patch("backend.layers.processing.upload_failures.app.trigger_slack_notification", mock_trigger_slack):
        event = {
            "dataset_version_id": "123",
            "collection_version_id": "456",
            "error": {},
            "execution_id": "arn:aws:states:us-west-2:123456789012:execution:MyStateMachine",
        }
        handle_failure(event, mock_context, delete_artifacts=False)
        mock_trigger_slack.assert_called_once()


def test_get_failure_slack_notification_message_with_dataset_version_id_none(
    get_collection_version_mock, sample_slack_header_block, sample_slack_status_block_empty, caplog
):
    dataset_version_id = None
    step_name = "Step 1"
    job_id = "123456"
    aws_regions = "us-west-2"
    execution_arn = "arn:aws:states:us-west-2:123456789012:execution:MyStateMachine"
    collection_version_id = "collection_version_id123"

    get_dataset_version_mock = Mock(return_value=None)
    get_business_logic_mock = Mock()
    get_business_logic_mock.get_dataset_version = get_dataset_version_mock
    get_business_logic_mock.get_collection_version = get_collection_version_mock
    get_business_logic_constructor_mock = Mock(return_value=get_business_logic_mock)

    with (
        patch(f"{module_path}.get_business_logic", get_business_logic_constructor_mock),
        caplog.at_level(logging.ERROR),
    ):
        result = get_failure_slack_notification_message(
            dataset_version_id, collection_version_id, step_name, job_id, aws_regions, execution_arn
        )
    assert result == {
        "blocks": [
            sample_slack_header_block,
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! Please follow the triage steps: https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit#bookmark=id.3ofm47y0709y\n"
                    "*Owner*: \n"
                    f"*Collection URL*: https://cellxgene.cziscience.com/collections/collection123\n"
                    f"*Collection Version URL*: https://cellxgene.cziscience.com/collections/{collection_version_id}\n"
                    "*Batch Job ID*: <https://us-west-2.console.aws.amazon.com/batch/v2/home?region=us-west-2"
                    "#jobs/detail/123456|123456>\n"
                    "*Step Function ARN*: "
                    "<https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/v2/executions"
                    "/details/arn:aws:states:us-west-2:123456789012:execution:MyStateMachine|arn:aws:states"
                    ":us-west-2:123456789012:execution:MyStateMachine>\n"
                    f"*Error Step*: {step_name}\n"
                    f"*Dataset Version ID*: None(not found)\n"
                    f"*Processing Status*:\n",
                },
            },
            sample_slack_status_block_empty,
        ]
    }
    assert "Dataset Version ID not found" in caplog.text


def test_get_failure_slack_notification_message_with_dataset_not_found(
    get_collection_version_mock, sample_slack_header_block, sample_slack_status_block_empty, caplog
):
    dataset_version_id = "dataset123"
    step_name = "Step 1"
    job_id = "123456"
    aws_regions = "us-west-2"
    execution_arn = "arn:aws:states:us-west-2:123456789012:execution:MyStateMachine"
    collection_version_id = "collection_version_id123"

    get_dataset_version_mock = Mock(return_value=None)
    get_business_logic_mock = Mock()
    get_business_logic_mock.get_dataset_version = get_dataset_version_mock
    get_business_logic_mock.get_collection_version = get_collection_version_mock
    get_business_logic_constructor_mock = Mock(return_value=get_business_logic_mock)

    with (
        patch(f"{module_path}.get_business_logic", get_business_logic_constructor_mock),
        caplog.at_level(logging.ERROR),
    ):
        result = get_failure_slack_notification_message(
            dataset_version_id, collection_version_id, step_name, job_id, aws_regions, execution_arn
        )

    assert result == {
        "blocks": [
            sample_slack_header_block,
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! Please follow the triage steps: https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit#bookmark=id.3ofm47y0709y\n"
                    "*Owner*: \n"
                    f"*Collection URL*: https://cellxgene.cziscience.com/collections/collection123\n"
                    f"*Collection Version URL*: https://cellxgene.cziscience.com/collections/{collection_version_id}\n"
                    "*Batch Job ID*: <https://us-west-2.console.aws.amazon.com/batch/v2/home?region=us-west-2"
                    "#jobs/detail/123456|123456>\n"
                    "*Step Function ARN*: "
                    "<https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/v2/executions"
                    "/details/arn:aws:states:us-west-2:123456789012:execution:MyStateMachine|arn:aws:states"
                    ":us-west-2:123456789012:execution:MyStateMachine>\n"
                    f"*Error Step*: {step_name}\n"
                    f"*Dataset Version ID*: {dataset_version_id}(not found)\n"
                    f"*Processing Status*:\n",
                },
            },
            sample_slack_status_block_empty,
        ]
    }
    assert "Dataset version ID dataset123 not found" in caplog.text
    get_dataset_version_mock.assert_called_with(DatasetVersionId(dataset_version_id))


def mock_collection_version(owner, version_id):
    MockCollectionVersion = Mock()
    MockCollectionVersion.owner = owner
    MockCollectionVersion.version_id = version_id
    return MockCollectionVersion


def test_get_failure_slack_notification_message_with_missing_collection(
    get_collection_version_mock, sample_slack_header_block, sample_slack_status_block, caplog
):
    dataset_version_id = "dataset123"
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
    get_business_logic_mock.get_collection_version = get_collection_version_mock
    get_business_logic_constructor_mock = Mock(return_value=get_business_logic_mock)

    with (
        patch(f"{module_path}.get_business_logic", get_business_logic_constructor_mock),
        caplog.at_level(logging.ERROR),
    ):
        result = get_failure_slack_notification_message(
            dataset_version_id, collection_version_id, step_name, job_id, aws_regions, execution_arn
        )

    assert result == {
        "blocks": [
            sample_slack_header_block,
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! Please follow the triage steps: https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit#bookmark=id.3ofm47y0709y\n"
                    f"*Owner*: \n"
                    f"*Collection URL*: https://cellxgene.cziscience.com/collections/{collection_id}\n"
                    f"*Collection Version URL*: https://cellxgene.cziscience.com/collections/{collection_version_id}\n"
                    "*Batch Job ID*: <https://us-west-2.console.aws.amazon.com/batch/v2/home?region=us-west-2"
                    "#jobs/detail/123456|123456>\n"
                    "*Step Function ARN*: "
                    "<https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/v2/executions"
                    "/details/arn:aws:states:us-west-2:123456789012:execution:MyStateMachine|arn:aws:states"
                    ":us-west-2:123456789012:execution:MyStateMachine>\n"
                    f"*Error Step*: {step_name}\n"
                    f"*Dataset Version ID*: {dataset_version_id}\n"
                    f"*Processing Status*:\n",
                },
            },
            sample_slack_status_block,
        ]
    }
    assert f"Collection {collection_id} not found" in caplog.text
    get_dataset_version_mock.assert_called_with(DatasetVersionId(dataset_version_id))
    get_unpublished_collection_version_from_canonical_mock.assert_called_with(CollectionId(collection_id))


def test_get_failure_slack_notification_message_with_dataset_and_collection(
    get_collection_version_mock, sample_slack_header_block, sample_slack_status_block
):
    dataset_version_id = "dataset123"
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
    get_business_logic_mock.get_collection_version = get_collection_version_mock
    get_business_logic_constructor_mock = Mock(return_value=get_business_logic_mock)

    with patch(f"{module_path}.get_business_logic", get_business_logic_constructor_mock):
        result = get_failure_slack_notification_message(
            dataset_version_id, collection_version_id, step_name, job_id, aws_regions, execution_arn
        )

    assert result == {
        "blocks": [
            sample_slack_header_block,
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": f"Dataset processing job failed! Please follow the triage steps: https://docs.google.com/document/d/1n5cngEIz-Lqk9737zz3makXGTMrEKT5kN4lsofXPRso/edit#bookmark=id.3ofm47y0709y\n"
                    f"*Owner*: {owner}\n"
                    f"*Collection URL*: https://cellxgene.cziscience.com/collections/{collection_id}\n"
                    f"*Collection Version URL*: https://cellxgene.cziscience.com/collections/{collection_version_id}\n"
                    "*Batch Job ID*: <https://us-west-2.console.aws.amazon.com/batch/v2/home?region=us-west-2"
                    "#jobs/detail/123456|123456>\n"
                    "*Step Function ARN*: "
                    "<https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/v2/executions"
                    "/details/arn:aws:states:us-west-2:123456789012:execution:MyStateMachine|arn:aws:states"
                    ":us-west-2:123456789012:execution:MyStateMachine>\n"
                    f"*Error Step*: {step_name}\n"
                    f"*Dataset Version ID*: {dataset_version_id}\n"
                    f"*Processing Status*:\n",
                },
            },
            sample_slack_status_block,
        ]
    }
    get_dataset_version_mock.assert_called_with(DatasetVersionId(dataset_version_id))
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
def dataset_version_id() -> str:
    return "example_dataset"


@pytest.fixture
def mock_delete_atac_fragment_files(monkeypatch) -> Mock:
    mock = Mock()
    monkeypatch.setattr(f"{module_path}.delete_atac_fragment_files", mock)
    return mock


@pytest.mark.usefixtures("mock_delete_atac_fragment_files")
class TestCleanupArtifacts:
    def test_cleanup_artifacts__validate_anndata_OK(
        self, mock_env_vars, mock_delete_many_from_s3, dataset_version_id, mock_delete_atac_fragment_files
    ):
        """Check that all artifacts are deleted for the given cases."""
        error_step = "validate_anndata"
        cleanup_artifacts(dataset_version_id, error_step)

        # Assertions
        mock_delete_atac_fragment_files.assert_not_called()
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["ARTIFACT_BUCKET"], dataset_version_id + "/")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["DATASETS_BUCKET"], dataset_version_id + ".")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["CELLXGENE_BUCKET"], dataset_version_id + ".cxg/")
        assert mock_delete_many_from_s3.call_count == 3

    def test_cleanup_artifacts__validate_atac_OK(
        self, mock_env_vars, mock_delete_many_from_s3, dataset_version_id, mock_delete_atac_fragment_files
    ):
        """Check that all artifacts are deleted for the given cases."""
        error_step = "validate_atac"
        cleanup_artifacts(dataset_version_id, error_step)

        # Assertions
        mock_delete_atac_fragment_files.assert_called_once_with(dataset_version_id)
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["DATASETS_BUCKET"], dataset_version_id + ".")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["CELLXGENE_BUCKET"], dataset_version_id + ".cxg/")
        assert mock_delete_many_from_s3.call_count == 2

    def test_cleanup_artifacts__None_OK(
        self, mock_env_vars, mock_delete_many_from_s3, dataset_version_id, mock_delete_atac_fragment_files
    ):
        """Check that all artifacts are deleted for the given cases."""
        error_step = None
        cleanup_artifacts(dataset_version_id, error_step)

        # Assertions
        mock_delete_atac_fragment_files.assert_called_once_with(dataset_version_id)
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["ARTIFACT_BUCKET"], dataset_version_id + "/")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["DATASETS_BUCKET"], dataset_version_id + ".")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["CELLXGENE_BUCKET"], dataset_version_id + ".cxg/")
        assert mock_delete_many_from_s3.call_count == 3

    def test_cleanup_artifacts__not_download_validate(
        self, mock_env_vars, mock_delete_many_from_s3, dataset_version_id, *args
    ):
        """Check that file in the artifact bucket are not delete if error_step is not download-validate."""
        cleanup_artifacts(dataset_version_id, "not_download_validate")

        # Assertions
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["DATASETS_BUCKET"], dataset_version_id + ".")
        mock_delete_many_from_s3.assert_any_call(mock_env_vars["CELLXGENE_BUCKET"], dataset_version_id + ".cxg/")
        assert mock_delete_many_from_s3.call_count == 2

    @patch.dict(os.environ, clear=True)
    def test_cleanup_artifacts__no_buckets(self, caplog, mock_delete_many_from_s3, dataset_version_id, *args):
        """Check that no files are deleted if buckets are not specified."""
        cleanup_artifacts(dataset_version_id)

        # Assertions
        mock_delete_many_from_s3.assert_not_called()
        assert FAILED_ARTIFACT_CLEANUP_MESSAGE in caplog.text
        assert FAILED_CXG_CLEANUP_MESSAGE in caplog.text
        assert FAILED_DATASET_CLEANUP_MESSAGE in caplog.text

    def test_cleanup_artifacts__delete_many_from_s3_error(
        self, caplog, mock_env_vars, mock_delete_many_from_s3, dataset_version_id, *args
    ):
        """Check that delete_many_from_s3 errors are logged but do not raise exceptions."""
        mock_delete_many_from_s3.side_effect = Exception("Boom!")
        cleanup_artifacts(dataset_version_id)

        # Assertions
        assert FAILED_ARTIFACT_CLEANUP_MESSAGE in caplog.text
        assert FAILED_CXG_CLEANUP_MESSAGE in caplog.text
        assert FAILED_DATASET_CLEANUP_MESSAGE in caplog.text


@pytest.fixture
def mock_dataset_version():
    dv = Mock()
    dv.artifacts = [
        DatasetArtifact(id=DatasetArtifactId(), uri="s3://bucket/uri", type=DatasetArtifactType.ATAC_INDEX),
        DatasetArtifact(id=DatasetArtifactId(), uri="s3://bucket/uri", type=DatasetArtifactType.ATAC_FRAGMENT),
    ]
    dv.status.atac_status = DatasetConversionStatus.UPLOADED
    return dv


class TestDeleteAtacFragmentFiles:
    @pytest.mark.parametrize(
        "atac_status", [DatasetConversionStatus.COPIED, DatasetConversionStatus.SKIPPED, DatasetConversionStatus.NA]
    )
    def test_delete_skipped(self, mock_business_logic, atac_status, mock_dataset_version):
        # Arrange
        dataset_version_id = "example_dataset"
        mock_dataset_version.status.atac_status = atac_status
        mock_business_logic.get_dataset_version.return_value = mock_dataset_version

        # Act
        delete_atac_fragment_files(dataset_version_id)

        # Assert
        mock_business_logic.get_atac_fragment_uris_from_dataset_version.assert_not_called()

    def test_delete_atac_fragment_files__OK(
        self, mock_delete_many_from_s3, mock_env_vars, mock_business_logic, mock_dataset_version, dataset_version_id
    ):
        """Check that atac fragment files are deleted."""
        # Arrange
        mock_business_logic.get_dataset_version.return_value = mock_dataset_version
        test_uris = ["uri1", "uri2"]
        mock_business_logic.get_atac_fragment_uris_from_dataset_version.return_value = ["uri1", "uri2"]

        # Act
        delete_atac_fragment_files(dataset_version_id)

        # Assertions
        for uri in test_uris:
            mock_delete_many_from_s3.assert_any_call(
                mock_env_vars["DATASETS_BUCKET"], os.path.join(os.environ.get("REMOTE_DEV_PREFIX", ""), uri)
            )

    def test_catch_errors(
        self,
        caplog,
        mock_delete_many_from_s3,
        mock_env_vars,
        mock_business_logic,
        mock_dataset_version,
        dataset_version_id,
    ):
        # Arrange
        mock_delete_many_from_s3.side_effect = Exception("Boom!")
        mock_business_logic.get_atac_fragment_uris_from_dataset_version.return_value = ["uri1"]
        # Act
        delete_atac_fragment_files(dataset_version_id)
        # Assert
        assert FAILED_ATAC_DATASET_MESSAGE[:-3] in caplog.text
