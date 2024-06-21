"""
This script is used to get the logs of a request from AWS cloudwatch given the request_id

It can also be used standal by modifying the constants in the __main__ block. This is useful for testing as well as
getting logs from an rdev environment which the cxg_admin CLI does not support.
"""

import json
import time
from typing import List

import boto3


def get(request_id: str, hours: int, deployment_stage: str, stackname: str) -> List[dict]:
    """
    Get the requests from from AWS cloudwatch given the request_id
    :param request_id:
    :param deployment_stage:
    :param stackname:
    """
    request_id = request_id.strip()
    client = boto3.client("logs")
    filter_pattern = f'{{$.request_id="{request_id}"}}'
    log_group_name = f"/dp/{deployment_stage}/{stackname}/backend"
    end_time = int(time.time()) * 1000
    start_time = end_time - hours * 60 * 60 * 1000

    def get_logs(next_token: str = None):
        kwargs = dict(
            logGroupName=log_group_name,
            logStreamNamePrefix="fargate/web/",
            filterPattern=filter_pattern,
            limit=100,
            unmask=False,
            startTime=start_time,
            endTime=end_time,
        )
        if next_token:
            kwargs["nextToken"] = next_token
        return client.filter_log_events(**kwargs)

    logs = get_logs()
    result = []
    beginning = False
    end = False
    while not all([beginning, end]):
        if logs["events"]:
            for log in logs["events"]:
                message = json.loads(log["message"])
                if message.get("type") == "RESPONSE":
                    end = True
                if message.get("type") == "REQUEST":
                    beginning = True
                result.append(message)
        if "nextToken" not in logs:
            break
        logs = get_logs(logs["nextToken"])
    return result


if __name__ == "__main__":
    import os

    os.environ["AWS_PROFILE"] = "single-cell-prod"
    request_id = "1407b9b9-c202-4b71-adbb-b08193379e3f"
    deployment_stage = "prod"
    stack_name = "prodstack"
    print(json.dumps(get(request_id, 3, deployment_stage, stack_name), indent=4))
