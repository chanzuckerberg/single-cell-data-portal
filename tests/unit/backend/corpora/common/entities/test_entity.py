import unittest

from backend.corpora.common.corpora_orm import DbDatasetProcessingStatus, UploadStatus
from backend.corpora.common.entities.dataset import Entity


class TestEntity(unittest.TestCase):
    def test__create_sub_object(self):
        test_params = {"row": {"upload_status": UploadStatus.WAITING}, "db_table": DbDatasetProcessingStatus}

        with self.subTest(test_params):
            result = Entity._create_sub_object(**test_params)
            self.assertIsInstance(result, DbDatasetProcessingStatus)
            self.assertEqual(UploadStatus.WAITING, result.upload_status)

        test_params.update(add_columns=dict(dataset_id="test_dataset_id"))
        with self.subTest(test_params):
            result = Entity._create_sub_object(**test_params)
            self.assertIsInstance(result, DbDatasetProcessingStatus)
            self.assertEqual("test_dataset_id", result.dataset_id)
            self.assertEqual(UploadStatus.WAITING, result.upload_status)

        test_params = {"row": {"fake_row": UploadStatus.WAITING}, "db_table": DbDatasetProcessingStatus}
        with self.subTest(test_params):
            with self.assertRaises(TypeError):
                Entity._create_sub_object(**test_params)
