from typing import List,Any, Union, Dict
from datetime import datetime, timezone

from config import get_sfn_client
from utils import utc_to_dt
from sfn import get_sfn_arn

# get execution results by step-function-ARN
def get_sfn_executions_by_arn(sfn_arn, status:str, start=None, end=None) -> List[Any]:
    client = get_sfn_client()
    executions = []
    has_next = True
    r = {} # mock an empty response
    
    # get ALL of the executions, as the API is limited to 1,000 per call.
    while has_next:
        kwargs = {"stateMachineArn":sfn_arn, "statusFilter":status, "maxResults":1_000}
        if r.get("nextToken",None):
            kwargs = kwargs | {"nextToken":r["nextToken"]}
        r = client.list_executions(**kwargs)
        executions.extend(r["executions"])
        has_next = r.get("nextToken",False)

    start = start if start else datetime(2025, 1, 1, 0, 0, 0)
    end = end if end else datetime(2025, 12, 31, 0, 0, 0)
    if isinstance(start,int):
        start = utc_to_dt(start)
    if isinstance(end, int):
        end = utc_to_dt(end)

    # NOTE:executions return datetime objects (in boto3)
    # filter executions by datetime range
    return [e for e in executions if e["startDate"] >= start and e["stopDate"] <= end]


# get the executions by name, single status only
def get_sfn_executions_by_name(sfn_name, status:str, start=None, end=None) -> List[Any]:
    arn = get_sfn_arn(sfn_name)
    return get_sfn_executions_by_arn(arn, status, start, end)


# get the executions by name, multiple status values
def get_all_sfn_executions_by_name(sfn_name, start=None, end=None, status:List[str]=["SUCCEEDED","FAILED"]) -> List[Any]:
    executions = []
    arn = get_sfn_arn(sfn_name)
    for _status in  status:
        executions.extend(get_sfn_executions_by_arn(arn, _status, start, end))
    return executions

# get the history of a specific execution
def get_sfn_execution_history(execution: Union[str,Dict[str,str]]) -> List[Any]:
    client = get_sfn_client()
    arn = execution if isinstance(execution,str) else execution.get("executionArn",None)
    return client.get_execution_history(executionArn=arn,reverseOrder=False)["events"]