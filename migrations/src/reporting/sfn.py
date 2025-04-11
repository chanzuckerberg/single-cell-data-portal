from typing import Union, List
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import re

import pandas as pd
from plotly import express as px

from config import get_sfn_client
from sfn.executions import get_all_sfn_executions_by_name, get_sfn_execution_history
from utils import ms_to_hours

@dataclass
class JobExecution:
    startedAt: datetime
    stoppedAt: datetime
    DATASET_VERSION_ID: str
    CURRENT_DATASET_VERSION_ID: Union[str,None]
    STEP_NAME: str
    status: str


# extract the current_dataset_version_id based on information in the job/step environment variables
def _pluck_current_dataset_version_id(environment_variables) -> str:
    '''
    DATASET_VERSION_ID: target id for the new dataset revision
    We want the original (to be migrated/validated) dataset version id. We'll call this
    CURRENT_DATASET_VERSION_ID.

    We get this from the MANIFEST environment variable, passed to a step. From there the CURRENT_DATASET_VERSION_ID
    is the name of the directory in the ARTIFACT_BUCKET e.g. s3://corpora-data-prod/__CURRENT_DATASET_VERSION_ID__.
    '''
    manifest = environment_variables.get("MANIFEST",None)
    # if not present, ignore
    if not manifest:
        return None

    # the value is a string; make it a json object
    manifest = json.loads(manifest)

    # get the anndata
    anndata_location = manifest.get("anndata",None)
    if not anndata_location:
        return None

    # get the ARTIFACT_BUCKET
    artifact_bucket = environment_variables.get("ARTIFACT_BUCKET") # this must exist if we get this far
    submission_bucket = "cellxgene-dataset-submissions/super" # alternative location for validation steps
    _prefix = "|".join([artifact_bucket, submission_bucket])
    
    r = re.search(rf"({_prefix})/(?P<cdvid>[0-9A-z-]+)", anndata_location)
    return r.group("cdvid")
    

# create JobExecution object from a successful sfn event
def _from_successful_task(event) -> JobExecution:
    data = json.loads(event["taskSucceededEventDetails"]["output"])
    env = {e["Name"]:e["Value"] for e in data["Container"]["Environment"]}
    
    return JobExecution(
        startedAt=data.get("StartedAt", data["CreatedAt"]),
        stoppedAt=data["StoppedAt"],
        DATASET_VERSION_ID=env["DATASET_VERSION_ID"],
        CURRENT_DATASET_VERSION_ID= _pluck_current_dataset_version_id(env),
        STEP_NAME = env["STEP_NAME"],
        status=data["Status"]
    )


# create JobExecution object from a failed sfn event
def _from_failed_task(event) -> JobExecution:
    data = json.loads(event["taskFailedEventDetails"]["cause"])
    env = {e["Name"]:e["Value"] for e in data["Container"]["Environment"]}
    return JobExecution(
        startedAt=data.get("StartedAt", data["CreatedAt"]),
        stoppedAt=data["StoppedAt"],
        DATASET_VERSION_ID=env["DATASET_VERSION_ID"],
        CURRENT_DATASET_VERSION_ID= _pluck_current_dataset_version_id(env),
        STEP_NAME = env["STEP_NAME"],
        status=data["Status"]
    )


# extract JobExecution objects from the step function history (only Successfull and Failed events)
def process_sfn_history(history) -> List[JobExecution]:
    records = []
    for event in history:
        if event["type"] == "TaskFailed":
            records.append(_from_failed_task(event))
        elif event["type"]=="TaskSucceeded":
            records.append(_from_successful_task(event))
    return records
        

# extract JobExecution objects from all executions of a step-function in a given time range
def batch_process_sfn(sfn_name:str, start:datetime, end:datetime, status: List[str]=["SUCCEEDED","FAILED"]):
    '''
    First get all executions in the time range.
    Then process the execution-history of each execution.
    Return Job Execution Reports
    '''
    print("retrieving executions")
    executions = get_all_sfn_executions_by_name(sfn_name, start, end,status)
    print("--DONE")

    print("retrieving histories")
    N = len(executions)
    records = []
    for i,e in enumerate(executions):
        if i%50==0:
            print(f"index={i} of {N}")
        history = get_sfn_execution_history(e)
        records.extend(process_sfn_history(history))
    print("--DONE")
    return records


# backfill CURRENT_DATASET_VERSION_ID
def _backfill_dsvid(df:pd.DataFrame) -> pd.DataFrame:
    # map DATASET -> CURRENT_DATASET
    map_ds = {r["DATASET_VERSION_ID"]:r["CURRENT_DATASET_VERSION_ID"] for _,r in df.iterrows() if r["CURRENT_DATASET_VERSION_ID"] is not None}

    # create new column of CURRENT_DATASET_VERSION_ID based on the above mapping
    CDVID = [map_ds.get(dvid,None) for dvid in df["DATASET_VERSION_ID"].values]

    df["CURRENT_DATASET_VERSION_ID"] = CDVID
    return df


# backfill durations
def _backfill_duration(df:pd.DataFrame) -> pd.DataFrame:
    #startedAt_dev	stoppedAt_dev
    _start = "startedAt"
    _stop = "stoppedAt"

    # row-wise duration (in ms)
    def duration(r):
        return r[_stop] - r[_start]
    
    df = df.assign(duration=df.apply(duration,axis=1))
    return df


# backfill various columns
def backfill(df:pd.DataFrame) -> pd.DataFrame:
    df = _backfill_dsvid(df)
    df = _backfill_duration(df)
    return df


# get dataframe with the summary for every step function execution
def describe_sfn_executions(sfn_name:str, start:datetime, end:datetime, status: List[str]=["SUCCEEDED","FAILED"]) -> pd.DataFrame:
    '''
    Get all executions of a step-function in a given time range.
    '''
    # get details
    records = batch_process_sfn(sfn_name, start, end, status)
    df = pd.DataFrame.from_records([r.__dict__ for r in records])
    df = backfill(df)
    return df

#---------------------------- COMPARISONS ----------------------------
'''
Reporting on differences between two environments (or other groups). 
'''

# row-wise computation of differences in duration (right - left)
def _compute_delta_duration(r, left_label:str, right_label:str) -> int:
    # ORDER is important
    return r[f"duration{right_label}"]- r[f"duration{left_label}"]


# backfill duration differences for every execution in the execution table
def _backfill_delta_duration(df:pd.DataFrame, left_label: str, right_label: str) -> pd.DataFrame:
    df = df.assign(delta_duration=df.apply(lambda x: _compute_delta_duration(x, left_label, right_label), axis=1))
    return df


# join on CURRENT_DATASET_VERSION_ID
def join_on_cdvid(left:pd.DataFrame, right:pd.DataFrame, left_label:str, right_label:str) -> pd.DataFrame:
    join_column = 'CURRENT_DATASET_VERSION_ID'
    return pd.merge(left,right, how="inner",on=join_column, suffixes=[f"_{left_label}",f"_{right_label}"])


# backfill joined dataset
def backfill_joined(df:pd.DataFrame, left_label: str, right_label: str) -> pd.DataFrame:
    df = _backfill_delta_duration(df, left_label, right_label)
    return df

# create comparison dataset
def create_comparison_dataset(left:pd.DataFrame, right:pd.DataFrame, left_label:str, right_label:str) -> pd.DataFrame:
    # join together into one table
    df = join_on_cdvid(left, right, left_label, right_label)

    # backfill the joined dataset with useful columns for comparison
    df = backfill_joined(df, left_label, right_label)


# COMPARISON: plot waterfall of duration differences for a single job (in two different environments)
def plot_job_duration_differences(left: pd.DataFrame, right: pd.DataFrame, left_label: str, right_label:str, job_name:str, count:int=50) -> px.Figure:
    '''
    Waterfall plot of duration differences.
    df: must be pre-joined on CURRENT_DATASET_VERSION_ID
    '''
    # create the comparison dataset
    df = create_comparison_dataset(left, right, left_label, right_label)

    # filter by job name
    df = df[df["STEP_NAME"]==job_name]

    # sort values in descending order, take top N
    df = df.sort_values("delta_duration",key=abs,ascending=False).iloc[0:count]

    # return the plot
    return px.bar(
        df,
        x="CURRENT_DATASET_VERSION_ID",
        y=df["delta_duration"].map(ms_to_hours),
        color=df["delta_duration"].map(lambda x: right_label if x>0 else left_label),
        title=f"Duration differences for {job_name}",
        labels={"delta_duration":"Duration difference (ms)"},
        text="delta_duration",
        orientation="v",   
    )