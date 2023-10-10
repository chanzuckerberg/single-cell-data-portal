import logging
from typing import Dict, List, Tuple
from unittest.mock import Mock

import pytest

from backend.layers.common.entities import DatasetProcessingStatus
from backend.layers.processing.publish_revisions import PublishRevisions


@pytest.fixture
def publish_revisions() -> PublishRevisions:
    business_logic = Mock()
    Mock()
    publish_revisions = PublishRevisions(business_logic)
    publish_revisions.schema_version = "1.0.0"
    return publish_revisions


@pytest.fixture
def publish_revisions_and_collections(
    tmpdir, publish_revisions, published_collection, revision, private
) -> Tuple[PublishRevisions, Dict[str, List]]:
    db = {
        published_collection.version_id.id: published_collection,
        revision[0].version_id.id: revision[0],
        revision[1].version_id.id: revision[1],
        private.version_id.id: private,
    }

    def _get_collection_version(collection_version_id):
        return db.get(collection_version_id.id)

    def _get_collections(filter):
        return db.values()

    publish_revisions.business_logic.get_collection_version = _get_collection_version
    publish_revisions.business_logic.get_collections = _get_collections
    return publish_revisions, {"published": [published_collection], "revision": revision, "private": [private]}


class TestPublishRevisions:
    def test_cleanup(self, publish_revisions, caplog):
        """Positive test case for cleanup method."""
        # Mock
        caplog.set_level(logging.INFO)
        publish_revisions.object_keys_to_delete = ["key1", "key2"]

        # Call
        publish_revisions.cleanup()

        # Assert
        assert any("Deleting migrated.h5ad files from s3." in rec.message for rec in caplog.records)
        publish_revisions.s3_provider.delete_files.assert_called_once_with(
            publish_revisions.artifact_bucket, ["key1", "key2"]
        )

    def test_check_datasets_no_errors(self, publish_revisions_and_collections):
        """positive test case for check_datasets method."""

        # Mock
        publish_revisions, collections = publish_revisions_and_collections
        published, revision = collections["revision"][0], collections["revision"][1]
        publish_revisions.business_logic.get_published_collection_version.return_value = published

        # Call
        errors = publish_revisions.check_datasets(revision)

        # Assert
        assert not errors
        assert publish_revisions.object_keys_to_delete == [
            f"{dataset.version_id.id}/migrated.h5ad" for dataset in published.datasets
        ]

    def test_check_datasets_not_latest_version(self, publish_revisions_and_collections):
        """Check datasets when one of the dataset has not been migrated to the latest schema version."""
        # Mock
        publish_revisions, collections = publish_revisions_and_collections
        published, revision = collections["revision"][0], collections["revision"][1]
        revision.datasets[0].metadata.schema_version = "0.0.1"
        publish_revisions.business_logic.get_published_collection_version.return_value = published

        # Call
        errors = publish_revisions.check_datasets(revision)

        # Assert
        assert errors[0]["message"] == "Dataset is not the latest schema version."
        assert publish_revisions.object_keys_to_delete == [f"{published.datasets[1].version_id.id}/migrated.h5ad"]

    def test_check_datasets_fail_status(self, publish_revisions_and_collections):
        """Check datasets when one of the dataset has failed during the ingestion process."""
        # Mock
        publish_revisions, collections = publish_revisions_and_collections
        published, revision = collections["revision"][0], collections["revision"][1]
        revision.datasets[0].status.processing_status = DatasetProcessingStatus.FAILURE
        error_message = "Error Message"
        revision.datasets[0].status.validation_message = error_message
        publish_revisions.business_logic.get_published_collection_version.return_value = published

        # Call
        errors = publish_revisions.check_datasets(revision)

        # Assert
        assert errors[0]["message"] == error_message
        assert publish_revisions.object_keys_to_delete == [f"{published.datasets[1].version_id.id}/migrated.h5ad"]

    def test_get_published_dataset_ids_from_collection_id(self, publish_revisions):
        """Positive test case for get_published_dataset_ids_from_collection_id method."""
        # Mock
        mock_collection_version = Mock()
        mock_dataset1 = Mock()
        mock_dataset1.dataset_id.id = "dataset1_id"
        mock_dataset1.version_id.id = "version1_id"
        mock_dataset2 = Mock()
        mock_dataset2.dataset_id.id = "dataset2_id"
        mock_dataset2.version_id.id = "version2_id"
        mock_collection_version.datasets = [mock_dataset1, mock_dataset2]
        publish_revisions.business_logic.get_published_collection_version.return_value = mock_collection_version

        # Call
        dataset_versions = publish_revisions.get_published_dataset_ids_from_collection_id(mock_collection_version)

        # Assert
        assert dataset_versions == {"dataset1_id": "version1_id", "dataset2_id": "version2_id"}

    def test_get_published_dataset_ids_from_collection_id__None(self, publish_revisions):
        """Test case for get_published_dataset_ids_from_collection_id method when the published version of the
        collection does not exists. This scenario should only be possible in test."""
        # Mock
        mock_collection_version = Mock()
        publish_revisions.business_logic.get_published_collection_version.return_value = None

        # Call
        dataset_versions = publish_revisions.get_published_dataset_ids_from_collection_id(mock_collection_version)

        # Assert the results
        assert dataset_versions == {}

    def test_run_pos(self, publish_revisions_and_collections, caplog):
        """Positive test case for run method."""
        # Mock necessary objects and methods
        caplog.set_level(logging.INFO)
        publish_revisions, collections = publish_revisions_and_collections
        _, revision = collections["revision"]
        publish_revisions.business_logic.get_collections.return_value = [collections]
        publish_revisions.check_datasets = Mock(return_value=[])

        # Call the method
        publish_revisions.run()

        # Assert the calls and behavior
        publish_revisions.check_datasets.assert_called_once_with(revision)
        assert "Publishing collection version." in caplog.messages
        publish_revisions.business_logic.publish_collection_version.assert_called_once_with(revision.version_id)
        assert "Deleting migrated.h5ad files from s3." in caplog.messages

    def test_run_neg(self, publish_revisions_and_collections, caplog):
        """Run method when check_datasets returns errors."""
        # Mock necessary objects and methods
        caplog.set_level(logging.INFO)
        publish_revisions, collections = publish_revisions_and_collections
        _, revision = collections["revision"]
        publish_revisions.business_logic.get_collections.return_value = [collections]
        publish_revisions.check_datasets = Mock(return_value=["error message"])

        # Call the method
        publish_revisions.run()

        # Assert the calls and behavior
        publish_revisions.check_datasets.assert_called_once_with(revision)
        assert "Unable to publish collection version." in caplog.messages
        publish_revisions.business_logic.publish_collection_version.assert_not_called()
        assert "Deleting migrated.h5ad files from s3." in caplog.messages


if __name__ == "__main__":
    pytest.main()
