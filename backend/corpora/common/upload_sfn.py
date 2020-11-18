import boto3
import json
import time

from .corpora_config import CorporaConfig

client = boto3.client("stepfunctions")


def start_upload_sfn(collection_uuid, dataset_uuid, url):
    input_parameters = {"collection_uuid": collection_uuid, "url": url, "dataset_uuid": dataset_uuid}
    sfn_name = f"{dataset_uuid}_{int(time.time())}"
    response = client.start_execution(
        stateMachineArn=CorporaConfig().upload_sfn_arn,
        # name=sfn_name,
        input=json.dumps(input_parameters),
    )
    return response
