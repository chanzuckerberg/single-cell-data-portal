from unittest import TestCase
from unittest.mock import patch, Mock

from backend.corpora.common.corpora_orm import generate_uuid
from backend.corpora.common.utils.corpora_constants import CorporaConstants
from backend.corpora.common.utils.exceptions import CorporaException
from backend.corpora.dataset_submissions.app import dataset_submissions_handler


class TestDatasetSubmissions(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.collection_uuid = generate_uuid()
        cls.incoming_curator_tag = "my_dataset.h5ad"
        cls.user_name = "user_name"
        cls.dataset_uuid = "12341234-1234-1234-1234-123412341234"

    def _test_missing_fields(self, **kwargs):
        s3_event = create_s3_event(**kwargs)
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    def test__missing_curator_tag__raises_error(self):
        self._test_missing_fields(key=f"{self.user_name}/{self.collection_uuid}/")

    def test__missing_collection_uuid__raises_error(self):
        self._test_missing_fields(key=f"{self.user_name}/should_have_been_a_uuid/some_key")

    def test__missing_username__raises_error(self):
        self._test_missing_fields(key=f"{self.collection_uuid}/some_key")

    @patch("backend.corpora.dataset_submissions.app.get_dataset_info")
    def test__non_existent_collection__raises_error(self, mock_get_dataset_info):
        s3_event = create_s3_event(key=f"{self.user_name}/{self.collection_uuid}/{self.incoming_curator_tag}")
        mock_get_dataset_info.return_value = None, None
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    @patch("backend.corpora.dataset_submissions.app.get_dataset_info")
    def test__non_owner__raises_error(self, mock_get_dataset_info):
        s3_event = create_s3_event(key=f"{self.user_name}/{self.collection_uuid}/{self.incoming_curator_tag}")
        mock_get_dataset_info.return_value = "not_owner", self.dataset_uuid
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    @patch("backend.corpora.dataset_submissions.app.get_dataset_info")
    @patch("backend.corpora.dataset_submissions.app.upload")
    def test__super_curator__upload(self, mock_upload: Mock, mock_get_dataset_info: Mock):
        s3_event = create_s3_event(
            key=f"{CorporaConstants.SUPER_CURATOR_NAME}/{self.collection_uuid}/{self.incoming_curator_tag}"
        )
        mock_upload.return_value = None
        mock_get_dataset_info.return_value = self.user_name, self.dataset_uuid
        dataset_submissions_handler(s3_event, None)
        mock_upload.assert_called()

    @patch("backend.corpora.dataset_submissions.app.get_dataset_info")
    @patch("backend.corpora.dataset_submissions.app.upload")
    def test__owner__upload(self, mock_upload: Mock, mock_get_dataset_info: Mock):
        s3_event = create_s3_event(key=f"{self.user_name}/{self.collection_uuid}/{self.incoming_curator_tag}")
        mock_upload.return_value = None
        mock_get_dataset_info.return_value = self.user_name, self.dataset_uuid
        dataset_submissions_handler(s3_event, None)
        mock_upload.assert_called()


def create_s3_event(bucket_name: str = "some_bucket", key: str = "", size: int = 0) -> dict:
    """
    Returns an S3 event dictionary with only the keys that the dataset submissions handler cares about
    :param bucket_name:
    :param key:
    :param size:
    :return:
    """
    return {"Records": [{"s3": {"bucket": {"name": bucket_name}, "object": {"key": key, "size": size}}}]}
