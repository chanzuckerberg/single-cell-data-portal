'''
Querying Batch Jobs
'''
from typing import List, Any
from enum import Enum
from datetime import datetime, timezone
from config import get_batch_client, get_job_queue
from utils import dt_to_utc, utc_to_dt

class STATUS_ENUM(Enum):
    succeeded = "SUCCEEDED"
    failed = "FAILED"
    running = "RUNNING"
    pending = "PENDING"
    submitted = "SUBMITTED"
    starting = "STARTING"



_BATCH_JOB_STATUS_LIST = ["SUCCEEDED", "FAILED", "RUNNING", "PENDING", "SUBMITTED", "STARTING"]

# get all jobs that beging inside an optional time range
def get_all_jobs_of_status(queue:str, status:str,start=None, end=None):
    all_jobs = []
    has_next = True
    r = {}
    client = get_batch_client()
    while has_next:
        if r.get("nextToken",None):
            r = client.list_jobs(jobQueue=queue,jobStatus=status,nextToken=r.get("nextToken"))
        else:
            r = client.list_jobs(jobQueue=queue,jobStatus=status)
        all_jobs.extend(r["jobSummaryList"])
        has_next = r.get("nextToken",False)

    # bookends (utc ms)
    start = start if start else datetime(2025,1,1,0,0,0)
    end = end if end else datetime(2025,12,31,0,0,0)

    # format times as milliseconds UTC
    if isinstance(start, datetime):
        start = dt_to_utc(start)
    if isinstance(end, datetime):
        end = dt_to_utc(end)

    # Use 'createdAt' since startedAt won't catch jobs that don't start
    return [j for j in all_jobs if j["createdAt"]> start and j["createdAt"] < end]

# get all jobs of all status
def get_all_jobs(queue:str, start=None, end=None) -> List[Any]:
    all_jobs = []
    for status in _BATCH_JOB_STATUS_LIST:
        all_jobs.extend(get_all_jobs_of_status(queue, status, start, end))
    return all_jobs

