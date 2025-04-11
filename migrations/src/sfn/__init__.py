from typing import List, Any
from config import get_sfn_client


# get step-function ARN by name
def get_sfn_arn(sfn_name:str) -> str:
    client = get_sfn_client()
    matches = [m for m in client.list_state_machines()["stateMachines"] if m["name"]==sfn_name]
    if len(matches)==0:
        raise Exception(f"No state machine found with name:{sfn_name}")
    else:
        return matches[0]["stateMachineArn"]
 