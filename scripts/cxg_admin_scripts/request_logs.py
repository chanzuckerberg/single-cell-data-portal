import json
from typing import List

import boto3


def get(request_id: str, deployment_stage: str, stackname: str) -> List[dict]:
    """
    Get the requests from from AWS cloudwatch given the request_id
    :param request_id:
    :param deployment_stage:
    :param stackname:
    """

    client = boto3.client("logs")
    filter_pattern = f'{{$.request_id="{request_id}"}}'
    log_group_name = f"/dp/{deployment_stage}/{stackname}/backend"

    def get_logs(next_token: str = None):
        kwargs = dict(
            logGroupName=log_group_name,
            logStreamNamePrefix="fargate/web/",
            filterPattern=filter_pattern,
            limit=100,
            unmask=False,
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
        logs = get_logs(logs["nextToken"])
    return result


if __name__ == "__main__":
    import os

    os.environ["AWS_PROFILE"] = "single-cell-dev"
    request_id = "cf017e14-7a04-426f-91a7-bf81318c27a6"
    deployment_stage = "rdev"
    stack_name = "pr-7217"
    print(json.dumps(get(request_id, deployment_stage, stack_name), indent=4))
