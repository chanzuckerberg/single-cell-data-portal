'''
Utilities to monitor the progress and errors for cxg datasets.
This operates on the level of datasets and not jobs.
Jobs can be completed (e.g. dataset_migrate) but the dataset is still being processed (e.g. validate_anndata).

NOTE: Not all jobs are specific to datasets/collections. For example, 'gather_collections' spans the entire database for collections. These jobs are ignored in this workflow 
'''
from typing import Union, List
from datetime import datetime, timezone, timedelta

import pandas as pd
from numpy import nan
from plotly import express as px
from utils import dt_to_utc
from jobs.jobs import STATUS_ENUM

_JOBNAME_COLUMN = "jobName"
_STATUS_COL = "status"
_CREATED_COL = "createdAt"
_START_COL = "startedAt"
_STOP_COL = "stoppedAt"
_STEP_NAME_COL = "STEP_NAME"

_DATASET_DURATION_COL = "duration"
_ACTIVE_STEP_COL = "active_step_name"
_ACTIVE_STEP_DURATION_COL = "active_step_duration"
_LAST_STEP_COL = "last_step_name"
_LAST_STEP_DURATION_COL = "last_step_duration"
_LAST_STEP_REASON_COL = "last_step_reason"
_DATASET_JOB_NAMES = ['dataset_migration','validate_anndata','validate_atac','add_labels','cxg'] # ordered

# get a start time based on createAt or startedAt
def _get_start_time(df: pd.DataFrame) -> Union[int,None]:
    """
    This is to help with edge case reports from AWS where a job may be stopped wihout ever starting.
    Then, we'd want to use the createdAt time.
    """
    _createdAt = df[_CREATED_COL].min()
    _startedAt = df[_START_COL].min()

    if pd.isna(_startedAt) and pd.isna(_createdAt):
        # no usable start time
        return None
    elif pd.isna(_startedAt) and not pd.isna(_createdAt):
        # use craeted at - sometimes AWS will stop a job that never started
        return _createdAt
    else:
        return _startedAt
    
# select the jobs that are specific to singular datasets
def _filter_dataset_jobs(df:pd.DataFrame) -> pd.DataFrame:
    return df[df[_JOBNAME_COLUMN].isin(_DATASET_JOB_NAMES)]

# aggregate operator: collapse the relevant status
def _aggregate_status(df:pd.DataFrame) -> str:
    # get the unique status results
    group_status = set(df[_STATUS_COL].unique())

    # return the aggregate status
    if STATUS_ENUM.failed.value in group_status:
        # if one job failed, the dataset failed
        return STATUS_ENUM.failed.value
    elif STATUS_ENUM.running.value in group_status:
        # if one job is running (and not failed), the dataset is running
        return STATUS_ENUM.running.value
    elif len(group_status.difference({STATUS_ENUM.succeeded.value}))==0:
        # if all jobs succeed, the datata succeeded
        return STATUS_ENUM.succeeded.value
    else:
        # otherwise it's running
        return STATUS_ENUM.running.value

# aggregate operator: earliest start time
def _aggregate_start_time(df:pd.DataFrame) -> Union[int,None]:
    # get the earliest start time
    return _get_start_time(df)

# aggregate operator: latest stop time
def _aggregate_stop_time(df:pd.DataFrame) -> Union[int, None]:
    """
    if any job is unfinished, the dataset is unfinished
    """
    # get status
    status = _aggregate_status(df)

    if status in {STATUS_ENUM.succeeded.value, STATUS_ENUM.failed.value}:
        # dataset is finished, use most recently stopped job
        return df['stoppedAt'].max()
    else:
        return None

# aggregate operator: active step
def _aggregate_active_step(df:pd.DataFrame) -> List[str]:
    # get running steps
    running = df[df[_STATUS_COL] == STATUS_ENUM.running.value]
    return list(running[_STEP_NAME_COL].unique())

# aggregate operator: active step duration
def _aggregate_active_step_duration(df:pd.DataFrame) -> int:
    # get running steps
    dfr = df[df[_STATUS_COL] == STATUS_ENUM.running.value]

    # get the earliest start time
    start_time = dfr[_START_COL].min()

    # get the latest stop time
    if pd.isna(start_time):
        return 0
    else:
        now = dt_to_utc(datetime.now(timezone.utc))
        return int(now-start_time)

# aggregate operator: dataset duration
def _aggregate_duration(df:pd.DataFrame) -> int:
    # get the earliest start time
    start_time = _aggregate_start_time(df)
    stop_time = _aggregate_stop_time(df)

    if not start_time:
        # never started - return 0
        return 0
    elif not stop_time:
        # started, but hasn't stopped - return time to now
        return datetime.now(timezone.utc)-start_time
    else:
        # return the time difference
        return stop_time-start_time

# aggrgate operator: last step
def _aggregate_last_step(df:pd.DataFrame) -> str:
    # get status
    status = _aggregate_status(df)

    # determine last step (depending on status)
    if status in {STATUS_ENUM.failed.value}:
        # if failed, give the most recently started failed job
        return df[df[_STATUS_COL]==STATUS_ENUM.failed.value].sort_values(by=_START_COL, ascending=False).iloc[0][_STEP_NAME_COL]
    else:
        # give the most recently started successful job
        return df[df[_STATUS_COL]==STATUS_ENUM.succeeded.value].sort_values(by=_START_COL, ascending=False).iloc[0][_STEP_NAME_COL]

# aggregate operator: last step duration
def _aggregate_last_step_duration(df:pd.DataFrame) -> int:
    # the step name
    step = _aggregate_last_step(df)

    # the step row (as a df)
    step = df[df[_STEP_NAME_COL]==step]

    start_time = _get_start_time(step)
    stop_time = step[_STOP_COL].max()

    return stop_time - start_time

# aggregate operator: last step reason
def _aggregate_last_step_reason(df:pd.DataFrame) -> Union[str,None]:
    '''
    This is only useful for failed jobs
    '''
    # get status
    status = _aggregate_status(df)

    if status == STATUS_ENUM.failed.value:
        # get the most recently started failed job
        failed = df[df[_STATUS_COL]==STATUS_ENUM.failed.value].sort_values(by=_START_COL, ascending=False).iloc[0]

        # return the reason
        return failed['reason'] or failed["statusReason"]
    else:
        return None
    

# get only dataset-level jobs
def get_dataset_jobs(jobs:pd.DataFrame) -> pd.DataFrame:
    return _filter_dataset_jobs(jobs)

# create a dataframe to report dataset status
def report_datasets(df:pd.DataFrame) -> pd.DataFrame:
    # filter dataset jobs
    dfds = _filter_dataset_jobs(df)

    # aggrregate information at the dataset level
    records = []
    for DATASET_VERSION_ID, _dfdv in dfds.groupby('DATASET_VERSION_ID'):
        # let's make sure the dataset-id exists
        _ids = list(set([x for x in _dfdv['DATASET_ID'] if not pd.isna(x)]))

        # aggregate status
        records.append({
            'DATASET_ID': _ids[0] if len(_ids) else nan,
            'DATASET_VERSION_ID': DATASET_VERSION_ID,
            _STATUS_COL: _aggregate_status(_dfdv),
            _START_COL: _aggregate_start_time(_dfdv),
            _STOP_COL: _aggregate_stop_time(_dfdv),
            _DATASET_DURATION_COL: _aggregate_duration(_dfdv),
            _ACTIVE_STEP_COL: _aggregate_active_step(_dfdv),
            _ACTIVE_STEP_DURATION_COL: _aggregate_active_step_duration(_dfdv),
            _LAST_STEP_COL: _aggregate_last_step(_dfdv),
            _LAST_STEP_DURATION_COL: _aggregate_last_step_duration(_dfdv),
            _LAST_STEP_REASON_COL: _aggregate_last_step_reason(_dfdv)

        })

    return pd.DataFrame.from_records(records)

# report dataset success breakdown
def plot_dataset_status_counts(df: pd.DataFrame, show:bool = False):
    f = px.histogram(df, x="status",histfunc="count", log_y=True, title="dataset migration status")
    if show:
        f.show()
    return f