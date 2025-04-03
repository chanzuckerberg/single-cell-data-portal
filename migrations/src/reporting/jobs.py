from datetime import datetime
from typing import Dict

import pandas as pd
from plotly import express as px

from config import get_job_queue
from jobs.jobs import get_all_jobs, STATUS_ENUM
from jobs.details import get_structured_job_details


def ms_to_minutes(time):
    return time/(1000*60)

def ms_to_hours(time):
    return ms_to_minutes(time)/60

# return job data as a dataframe
def describe_migration_jobs(start_time:datetime, end_time: datetime) -> pd.DataFrame:
    # get details
    _QUEUE = get_job_queue()
    details = get_structured_job_details(get_all_jobs(_QUEUE, start_time, end_time))
    return pd.DataFrame.from_records(details)


# create data visualization of job counts by job name and status
def plot_job_status_counts(df:pd.DataFrame, show:bool=False):
    f = px.histogram(df, x="jobName",histfunc='count',color="status", title="status by job", log_y=True, barmode="group")
    if show:
        f.show()
    return f


# create data visualization of job counts by job name and status
def plot_job_durations(df:pd.DataFrame, show:bool=False):
    f = px.histogram(df,x="jobName",histfunc="avg",y=df["duration"].map(ms_to_hours), color="status",barmode='group', log_y=True, title="duration (hours) by job") 
    if show:
        f.show()
    return f
