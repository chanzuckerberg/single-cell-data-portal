from unittest.mock import Mock, patch

from backend.common.utils.exceptions import (
    CorporaException,
    NonExistentCollectionException,
    NonExistentDatasetException,
)
from backend.layers.common.entities import EntityId
from backend.layers.processing.submissions.app import dataset_submissions_handler
from tests.unit.backend.layers.common.base_test import BaseTest


class TestDatasetSubmissions(BaseTest):
    def setUp(self) -> None:
        super().setUp()
        self.user_name = "test_user_id"
        self.mock = patch(
            "backend.layers.processing.submissions.app.get_business_logic", return_value=self.business_logic
        )
        self.mock.start()

    def tearDown(self):
        self.mock.stop()

    def test__missing_curator_file_name__raises_error(self):
        mock_collection_id = EntityId()
        s3_event = create_s3_event(key=f"{self.user_name}/{mock_collection_id}/")
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    def test__missing_collection_id__raises_error(self):
        mock_collection_id = EntityId()
        mock_dataset_id = EntityId()

        s3_event = create_s3_event(key=f"{self.user_name}/{mock_collection_id}/{mock_dataset_id}")
        with self.assertRaises(NonExistentCollectionException):
            dataset_submissions_handler(s3_event, None)

    def test__nonexistent_dataset__raises_error(self):
        version = self.generate_unpublished_collection()
        mock_dataset_id = EntityId()

        s3_event = create_s3_event(key=f"{self.user_name}/{version.version_id}/{mock_dataset_id}")
        with self.assertRaises(NonExistentDatasetException):
            dataset_submissions_handler(s3_event, None)

    def test__non_owner__raises_error(self):
        version = self.generate_unpublished_collection(owner="someone_else")
        mock_dataset_id = EntityId()

        s3_event = create_s3_event(key=f"{self.user_name}/{version.version_id}/{mock_dataset_id}")
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    def test__upload_update_by_dataset_id_owner__OK(self):
        """
        Processing starts when an update of a dataset is uploaded by its ID by the collection owner.
        """
        version = self.generate_unpublished_collection()
        dataset_version_id, _ = self.business_logic.create_empty_dataset(version.version_id)

        mock_ingest = self.business_logic.ingest_dataset = Mock()

        s3_event = create_s3_event(key=f"{self.user_name}/{version.version_id}/{dataset_version_id}")
        dataset_submissions_handler(s3_event, None)
        mock_ingest.assert_called()

    def test__upload_update_by_dataset_id_super__OK(self):
        """
        Processing starts when an update of a dataset is uploaded by its ID by a super curator
        """
        version = self.generate_unpublished_collection()
        dataset_version_id, _ = self.business_logic.create_empty_dataset(version.version_id)

        mock_ingest = self.business_logic.ingest_dataset = Mock()

        s3_event = create_s3_event(key=f"super/{version.version_id}/{dataset_version_id}")
        dataset_submissions_handler(s3_event, None)
        mock_ingest.assert_called()

    def test__upload_update_by_dataset_canonical_id__OK(self):
        """
        Processing starts when an update of a dataset is uploaded using canonical ids

        """
        version = self.generate_unpublished_collection()
        _, dataset_id = self.business_logic.create_empty_dataset(version.version_id)

        mock_ingest = self.business_logic.ingest_dataset = Mock()

        s3_event = create_s3_event(key=f"{self.user_name}/{version.collection_id}/{dataset_id}")
        dataset_submissions_handler(s3_event, None)
        mock_ingest.assert_called()


def create_s3_event(bucket_name: str = "some_bucket", key: str = "", size: int = 0) -> dict:
    """
    Returns an S3 event dictionary with only the keys that the dataset submissions handler cares about
    :param bucket_name:
    :param key:
    :param size:
    :return:
    """
    return {"Records": [{"s3": {"bucket": {"name": bucket_name}, "object": {"key": key, "size": size}}}]}
