import unittest
from backend.corpora.common.upload_sfn import start_upload_sfn, client

class Test_Uploader_SFN(unittest.TestCase):
    def test_happy(self):
        response = start_upload_sfn("test_collection_id", "test_dataset_uuid", "test_url")
        client.describe_execution(response['executionArn'])
        self.assertTrue(True)
