import boto3
import json
import time

from .corpora_config import CorporaConfig
import os


_stepfunctions_client = None


def get_stepfunctions_client():
    global _stepfunctions_client
    if not _stepfunctions_client:
        _stepfunctions_client = boto3.client("stepfunctions", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"))
    return _stepfunctions_client


def start_upload_sfn(collection_uuid, dataset_uuid, url):
    input_parameters = {"collection_uuid": collection_uuid, "url": url, "dataset_uuid": dataset_uuid}
    sfn_name = f"{dataset_uuid}_{int(time.time())}"
    response = get_stepfunctions_client().start_execution(
        stateMachineArn=CorporaConfig().upload_sfn_arn,
        name=sfn_name,
        input=json.dumps(input_parameters),
    )
    return response
