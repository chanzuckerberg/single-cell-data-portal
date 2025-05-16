from typing import List
import json

import pandas as pd
from logs.logs import get_streams_logs

_LOG_RESULT_COLNAME = "logResults"
_LOG_FETCHED_COLNAME = "logFetched"

def filter_errors(events):
    return [e for e in events if e["levelname"]=="ERROR"]

def filter_info(events):
    return [e for e in events if e["levelname"]=="INFO"]

def extract_error(event):
    keys = ["message","lineno","pathname","exc_info"]
    return {k:event.get(k,None) for k in keys}

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



# get the logs for all jobs in the dataframe
def get_job_logs(df_jobs:pd.DataFrame):
    # initialize dataframe to hold log results
    _df_logs = pd.DataFrame()
    
    # group fetch queries by log-group
    for log_group, dfg in df_jobs.groupby("awslogs-group"):
        # pull streams
        log_streams = dfg['logStreamName'].values
    
        # get logs from AWS
        log_results = get_streams_logs(log_group, log_streams)
    
        # build the dataframe from dictionary
        _df = pd.DataFrame.from_dict({
            "logStreamName":[x[1] for x in log_results],
            _LOG_RESULT_COLNAME: [x[2] for x in log_results],
            _LOG_FETCHED_COLNAME: [x[3] for x in log_results]
        })

        # concatenate
        _df_logs = pd.concat([_df_logs, _df],axis=0)
    return _df_logs

# get the logs for all jobs in the dataframe
def merge_job_logs(df_jobs:pd.DataFrame) -> pd.DataFrame:
    # create a dataframe with new columns with data from logs
    if _LOG_RESULT_COLNAME not in df_jobs.columns:
        df_logs = get_job_logs(df_jobs)

        # return new dataframe with log-columns
        return df_jobs.merge(df_logs, on="logStreamName", how="left")
    else:
        return df_jobs

