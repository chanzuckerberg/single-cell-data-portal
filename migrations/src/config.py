from typing import Union
import os
from enum import Enum
from boto3 import Session

class STAGE(Enum):
    prod = 'prod'
    dev = 'dev'


_BATCH_CLIENT = None,  
_LOGS_CLIENT = None
_JOB_QUEUE = None

def _get_profile(stage:STAGE) -> str:
    return f"single-cell-{stage.value}"

# get the aws session object, using the appropriate profile
def _get_session(stage:STAGE) -> Session:
    # will error out if you don't have this profile set in your .aws/config
    return Session(profile_name=_get_profile(stage), region_name="us-west-2")


# set the environment stage against which to monitor jobs and logs
def set_stage(stage: Union[STAGE, str]=STAGE.dev) -> None:
    if isinstance(stage, str):
        try:
            stage = STAGE(stage)
        except:
            raise ValueError(f"The provided stage '{stage}' is not a valid stage.")
    
    # set the environment variable
    os.environ['MIGRATION_MONITOR_STAGE'] = stage.value

    # get aws session
    _session = _get_session(stage)

    # set global constants based on stage 
    global _BATCH_CLIENT, _LOGS_CLIENT, _JOB_QUEUE
    _BATCH_CLIENT = _session.client('batch')
    _LOGS_CLIENT = _session.client("logs")
    
    # set other constants
    _JOB_QUEUE = f"schema_migration-{stage.value}"

def get_batch_client():
    return _BATCH_CLIENT

def get_logs_client():
    return _LOGS_CLIENT

def get_job_queue():
    return _JOB_QUEUE