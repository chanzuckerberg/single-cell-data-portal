from unittest import TestCase
from unittest.mock import patch, Mock

from backend.common.corpora_orm import generate_id
from backend.common.utils.exceptions import (
    CorporaException,
    NonExistentCollectionException,
    NonExistentDatasetException,
)
from backend.dataset_pipeline.submissions.app import dataset_submissions_handler


class TestDatasetSubmissions(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.collection_id = generate_id()
        cls.user_name = "user_name"
        cls.dataset_id = "12341234-1234-1234-1234-123412341234"

    def _test_missing_fields(self, **kwargs):
        s3_event = create_s3_event(**kwargs)
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    def test__missing_curator_file_name__raises_error(self):
        self._test_missing_fields(key=f"{self.user_name}/{self.collection_id}/")

    @patch("backend.dataset_pipeline.submissions.app.Dataset.get")
    def test__missing_collection_id__raises_error(self, mock_dataset_get):
        mock_dataset_get.return_value = None
        s3_event = create_s3_event(key=f"{self.user_name}/{self.collection_id}/{self.dataset_id}")
        with self.assertRaises(NonExistentCollectionException):
            dataset_submissions_handler(s3_event, None)

    @patch("backend.dataset_pipeline.submissions.app.Collection.get")
    def test__nonexistent_dataset__raises_error(self, mock_collection_get):
        mock_collection_get.return_value = True
        s3_event = create_s3_event(key=f"{self.user_name}/{self.collection_id}/{self.dataset_id}")
        with self.assertRaises(NonExistentDatasetException):
            dataset_submissions_handler(s3_event, None)

    @patch("backend.dataset_pipeline.submissions.app.Dataset.get")
    def test__non_owner__raises_error(self, mock_dataset_get):
        mock_dataset_get.return_value = make_dataset_mock("now_owner", self.dataset_id)
        s3_event = create_s3_event(key=f"{self.user_name}/{self.collection_id}/{self.dataset_id}")
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    @patch("backend.dataset_pipeline.submissions.app.Dataset.get")
    def test__upload_new_by_dataset_id__Error(self, mock_get_dataset: Mock):
        """processing fails when a new dataset is uploaded using a dataset uuid that is not part of the collection."""
        uuid = self.dataset_id
        mock_get_dataset.return_value = None
        with self.assertRaises(CorporaException):
            s3_event = create_s3_event(key=f"{self.user_name}/{self.collection_id}/{uuid}")
            dataset_submissions_handler(s3_event, None)

    @patch("backend.dataset_pipeline.submissions.app.Dataset.get")
    def test__upload_update_by_dataset_id__OK(self, mock_get_dataset):
        """processing starts when an update of a dataset is uploaded by its ID."""
        mock_dataset = make_dataset_mock(self.user_name, self.dataset_id)
        mock_get_dataset.return_value = mock_dataset
        with self.subTest("owner"), patch("backend.dataset_pipeline.submissions.app.upload") as mock_upload:
            s3_event = create_s3_event(key=f"{self.user_name}/{self.collection_id}/{self.dataset_id}")
            dataset_submissions_handler(s3_event, None)
            mock_upload.assert_called()
        with self.subTest("super"), patch("backend.dataset_pipeline.submissions.app.upload") as mock_upload:
            s3_event = create_s3_event(key=f"super/{self.collection_id}/{self.dataset_id}")
            dataset_submissions_handler(s3_event, None)
            mock_upload.assert_called()


def make_dataset_mock(owner, dataset_id):
    mock_dataset = Mock()
    mock_dataset.collection.owner = owner
    mock_dataset.id = dataset_id
    return mock_dataset


def create_s3_event(bucket_name: str = "some_bucket", key: str = "", size: int = 0) -> dict:
    """
    Returns an S3 event dictionary with only the keys that the dataset submissions handler cares about
    :param bucket_name:
    :param key:
    :param size:
    :return:
    """
    return {"Records": [{"s3": {"bucket": {"name": bucket_name}, "object": {"key": key, "size": size}}}]}
