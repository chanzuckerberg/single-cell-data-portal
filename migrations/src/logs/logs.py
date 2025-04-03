from typing import List
import json
from config import get_logs_client

def _parse_log_events(events:List[str]):
    # we only care about logger messages
    messages = []
    for e in events:
        try:
            m = json.loads(e["message"])
        except:
            continue
        messages.append(m)
    return messages

# fetch logs for a single log-stream
def get_stream_logs(log_group: str, log_stream: str):
    client = get_logs_client()
    success = True
    _ptoken = '__'
    _ntoken = ''
    events = []
    count = 0
    max_count = 5
    while _ntoken != _ptoken and count < max_count:
        count+=1
        # keyword arguments
        kwargs = {"logGroupName":log_group, "logStreamName":log_stream, "startFromHead":True}
        kwargs = kwargs if _ntoken =='' else kwargs | {"nextToken":_ntoken}
        try:
            r = client.get_log_events(**kwargs)
            events.extend(r["events"])
            _ptoken = _ntoken
            _ntoken = r["nextForwardToken"]
        except Exception as e:
            success = False
            

    return (log_group, log_stream, _parse_log_events(events), success)


# fetch logs for multiple log-streams
def get_streams_logs(log_group: str, log_streams: List[str]):
    #TODO:[EM] this should be made async enabled with aiobotocore
    return [get_stream_logs(log_group, stream) for stream in log_streams]