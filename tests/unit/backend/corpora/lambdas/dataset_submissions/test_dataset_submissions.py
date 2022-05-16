from unittest import TestCase
from unittest.mock import patch

from backend.corpora.common.corpora_orm import generate_uuid
from backend.corpora.common.utils.exceptions import CorporaException
from backend.corpora.dataset_submissions.app import dataset_submissions_handler


class TestDatasetSubmissions(TestCase):
    def test__missing_curator_tag__raises_error(self):
        s3_event = create_s3_event(bucket_name="some_bucket", key=f"user_name/{generate_uuid()}/")
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    def test__missing_collection_uuid__raises_error(self):
        s3_event = create_s3_event(key="user_name/should_have_been_a_uuid/some_key")
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    def test__missing_username__raises_error(self):
        s3_event = create_s3_event(key="should_have_been_a_uuid/some_key")
        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)

    @patch("backend.corpora.dataset_submissions.app.get_dataset_info")
    def test__non_existent_collection__raises_error(self, mock_get_dataset_info):
        collection_uuid = generate_uuid()
        incoming_curator_tag = "my_dataset.h5ad"
        s3_event = create_s3_event(bucket_name="expected_bucket", key=f"{collection_uuid}/{incoming_curator_tag}")

        mock_get_dataset_info.return_value = None, None

        with self.assertRaises(CorporaException):
            dataset_submissions_handler(s3_event, None)


def create_s3_event(bucket_name: str = "some_bucket", key: str = "", size: int = 0) -> dict:
    """
    Returns an S3 event dictionary with only the keys that the dataset submissions handler cares about
    :param bucket_name:
    :param key:
    :param size:
    :return:
    """
    return {"Records": [{"s3": {"bucket": {"name": bucket_name}, "object": {"key": key, "size": size}}}]}
