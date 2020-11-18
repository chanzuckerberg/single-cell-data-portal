import unittest
import os
import json
from backend.corpora.common.upload_sfn import start_upload_sfn, client
from backend.corpora.common.corpora_config import CorporaConfig


class Test_Uploader_SFN(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        os.environ["UPLOAD_SFN_ARN"] = CorporaConfig(deployment="dev").upload_sfn_arn

    def test_happy(self):
        input_parameters = {
            "collection_uuid": "test_collection_id",
            "dataset_uuid": "test_dataset_uuid",
            "url": "test_url",
        }
        response = start_upload_sfn(**input_parameters)
        response = client.describe_execution(executionArn=response["executionArn"])
        self.assertDictEqual(input_parameters, json.loads(response["input"]))
