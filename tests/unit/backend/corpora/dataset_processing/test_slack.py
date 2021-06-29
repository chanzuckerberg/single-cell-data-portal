import json
import os

from backend.corpora.dataset_processing.slack import format_slack_message
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestDatasetProcessing(DataPortalTestCase):
    def test_format_slack_message(self):
        dataset = self.generate_dataset(self.session)
        os.environ["AWS_BATCH_JOB_ID"] = "test_job_id"
        message = format_slack_message(dataset.id)
        json.loads(message)
