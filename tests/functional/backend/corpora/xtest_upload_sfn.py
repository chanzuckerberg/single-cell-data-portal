import json
import os
import unittest

import boto3

from backend.corpora.common.upload_sfn import start_upload_sfn

if os.getenv("DEPLOYMENT_STAGE") == "test":
    step_definition = json.loads(
        """{
      "Comment": "A Hello World example demonstrating various state types of the Amazon States Language",
      "StartAt": "Pass",
      "States": {
        "Pass": {
          "Type": "Pass",
          "End": true
        }
      }
    }
    """
    )
    client_sfn = boto3.client("stepfunctions", endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None)
    client_ssm = boto3.client("secretsmanager", endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None)
    secret_id = "corpora/backend/test/config"
    old_secret = client_ssm.get_secret_value(SecretId=secret_id)["SecretString"]

    def setUpModule():
        response = client_sfn.create_state_machine(
            name="fake_step",
            definition=json.dumps(step_definition),
            roleArn="arn:aws:iam::123456789012:role/corpora-dataset-uploader-test-sfn-service",
        )
        secret_string = json.dumps({"upload_sfn_arn": response["stateMachineArn"]})
        client_ssm.update_secret(SecretId=secret_id, SecretString=secret_string)

    def tearDownModule():
        upload_sfn_arn = json.loads(client_ssm.get_secret_value(SecretId=secret_id)["SecretString"])["upload_sfn_arn"]
        client_sfn.delete_state_machine(stateMachineArn=upload_sfn_arn)
        client_ssm.update_secret(SecretId=secret_id, SecretString=old_secret)


class Test_Uploader_SFN(unittest.TestCase):
    def test__happy__OK(self):
        """
        Invokes the step functions and verified the correct parameters have been passed.
        :return:
        """
        client_sfn = boto3.client("stepfunctions", endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None)
        input_params = {"collection_uuid": "test_collection_id", "url": "test_url", "dataset_uuid": "test_dataset_uuid"}
        response = start_upload_sfn(**input_params)
        response = client_sfn.describe_execution(executionArn=response["executionArn"])
        self.assertDictEqual(input_params, json.loads(response["input"]))
