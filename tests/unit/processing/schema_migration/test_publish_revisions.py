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

    def test_run_pos__do_not_publish_list(self, publish_revisions_and_collections, caplog):
        """Run method when revision is skipped via do_not_publish_list."""
        # Mock necessary objects and methods
        caplog.set_level(logging.INFO)
        publish_revisions, collections = publish_revisions_and_collections
        _, revision = collections["revision"]
        publish_revisions.business_logic.get_collections.return_value = [collections]
        publish_revisions.do_not_publish_list = [revision.version_id.id]
        publish_revisions.check_datasets = Mock(return_value=[])

        # Call the method
        publish_revisions.run()

        # Assert the calls and behavior
        publish_revisions.check_datasets.assert_not_called()
        assert "Skipping collection version, it is in the do not publish list" in caplog.messages
        publish_revisions.business_logic.publish_collection_version.assert_not_called()


if __name__ == "__main__":
    pytest.main()
