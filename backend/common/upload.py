import json
import os
import time

import boto3

from backend.common.corpora_config import CorporaConfig

_stepfunctions_client = None


def get_stepfunctions_client():
    global _stepfunctions_client
    if not _stepfunctions_client:
        _stepfunctions_client = boto3.client("stepfunctions", endpoint_url=os.getenv("BOTO_ENDPOINT_URL") or None)
    return _stepfunctions_client


def start_upload_sfn(collection_id, dataset_id, url):
    input_parameters = {
        "collection_id": collection_id,
        "url": url,
        "dataset_id": dataset_id,
    }
    sfn_name = f"{dataset_id}_{int(time.time())}"
    response = get_stepfunctions_client().start_execution(
        stateMachineArn=CorporaConfig().upload_sfn_arn,
        name=sfn_name,
        input=json.dumps(input_parameters),
    )
    return response
