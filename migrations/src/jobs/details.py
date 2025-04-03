from datetime import datetime, timezone

import pandas as pd
from config import get_batch_client
from utils import dt_to_utc, utc_to_dt

# pluck job-description details as a flat dictionary
def _pluck_details(jobdescription):
    # top level keys
    keys = ['jobName','jobId','jobQueue','status','statusReason','createdAt','stoppedAt','startedAt','jobDefinition']
    d = {k:jobdescription.get(k,None) for k in keys}

    # attempt.container level keys (only pluck last attempt)
    attempt_container_keys = ["reason"]
    if len(jobdescription["attempts"]):
        attempt_container = jobdescription["attempts"][-1]["container"]
        d = d | {k:attempt_container.get(k,None) for k in attempt_container_keys}
    else:
        d = d | {k:pd.NA for k in attempt_container_keys}

    # container level keys
    container = jobdescription['container']
    container_keys = ['logStreamName']
    d = d | {k:container.get(k,None) for k in container_keys}

    # container.environment level keys
    env = container['environment']
    env_keys = ['DATASET_VERSION_ID','ARTIFACT_BUCKET','DATASET_ID','COLLECTION_VERSION_ID','COLLECTION_ID','DATASETS_BUCKET','STEP_NAME']
    d = d | {e['name']:e['value'] for e in env if e['name'] in env_keys}

    # container.logConfiguration.option keys
    logopts = container['logConfiguration']['options']
    logopts_keys = ['awslogs-group']
    d = d | {k:logopts.get(k,None) for k in logopts_keys}

    # compute durations
    if d['startedAt'] is not None and d['stoppedAt'] is not None:
        d['duration'] = d['stoppedAt'] - d['startedAt']
    elif d['startedAt'] is not None and d['stoppedAt'] is None:
        now = dt_to_utc(datetime.now(timezone.utc))
        d['duration'] = now - d['startedAt']
    else:
        d['duration'] = pd.NA

    return d


# get all job descriptions
def get_job_details(jobs):
    results = []
    print(f'Job Count: {len(jobs)}')
    client = get_batch_client()
    for i in range(0, len(jobs), 100):  # describe_jobs takes max 100
        chunk = [j["jobId"] for j in jobs[i:i+100]]
        try:
            response = client.describe_jobs(jobs=chunk)
        except Exception as e:
            print(e)
        results.extend(response["jobs"])
    return results


def get_structured_job_details(jobs):
    return [_pluck_details(d) for d in get_job_details(jobs)]