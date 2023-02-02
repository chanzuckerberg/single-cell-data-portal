import os

from backend.layers.processing.upload_failures.app import get_failure_slack_notification_message
from tests.unit.backend.layers.common.base_api_test import BaseAPIPortalTest


class TestDatasetProcessing(BaseAPIPortalTest):
    def test_format_slack_message(self):
        dataset = self.generate_dataset()
        os.environ["AWS_BATCH_JOB_ID"] = "test_job_id"
        message = get_failure_slack_notification_message(dataset.dataset_version_id)
        self.assertTrue(message["blocks"])
