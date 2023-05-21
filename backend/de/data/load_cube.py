import os
import shutil
import subprocess

import boto3

from backend.wmg.data.schemas.corpus_schema import (
    INTEGRATED_ARRAY_NAME,
    OBS_ARRAY_NAME,
    VAR_ARRAY_NAME,
)
from backend.wmg.data.utils import log_func_runtime

stack_name = os.environ.get("REMOTE_DEV_PREFIX")
wmg_bucket_name = os.environ.get("WMG_BUCKET")
artifact_bucket_name = os.environ.get("ARTIFACT_BUCKET")


def _get_wmg_bucket_path():
    if stack_name:
        return f"s3://{wmg_bucket_name}{stack_name}"
    else:
        return f"s3://{wmg_bucket_name}"


def make_snapshot_active(snapshot_id):
    """
    Update snapshot pointer and remove older datasets
    """
    write_snapshot_id(snapshot_id)
    remove_oldest_datasets(snapshot_id)


def remove_oldest_datasets(timestamp):
    """
    Remove all data stored under the oldest timestamp
    """
    s3 = boto3.resource("s3")
    wmg_bucket = s3.Bucket(wmg_bucket_name)
    objects = wmg_bucket.objects.filter(Prefix=stack_name.strip("/")) if stack_name else wmg_bucket.objects.all()

    candidate_to_delete = []
    for obj in objects:
        try:
            tokens = obj.key.split("/")
            timestamp_token = tokens[1] if stack_name else tokens[0]
            is_timestamp = timestamp_token[:10].isdigit()
            if is_timestamp:
                candidate_to_delete.append((timestamp_token, obj))
        except Exception:
            pass

    timestamps = sorted({x[0] for x in candidate_to_delete})

    if len(timestamps) > 2:
        timestamps_to_delete = list(timestamps)[:-2]
    else:
        return

    for timestamp, object in candidate_to_delete:
        if timestamp in timestamps_to_delete:
            object.delete()


@log_func_runtime
def upload_artifacts_to_s3(snapshot_path, s3_key):
    # cleanup files we do not wish to upload
    shutil.rmtree(f"{snapshot_path}/{INTEGRATED_ARRAY_NAME}")
    shutil.rmtree(f"{snapshot_path}/{OBS_ARRAY_NAME}")
    shutil.rmtree(f"{snapshot_path}/{VAR_ARRAY_NAME}")

    sync_command = ["aws", "s3", "sync", snapshot_path, f"{_get_wmg_bucket_path()}/{s3_key}"]
    subprocess.run(sync_command)


def write_snapshot_id(timestamp):
    """
    Update static timestamp in s3 to match the latest
    """
    s3 = boto3.resource("s3")
    file_key = "latest_snapshot_identifier"
    file_path = f"{stack_name.strip('/')}/{file_key}" if stack_name else f"{file_key}"
    object = s3.Object(wmg_bucket_name, file_path)
    object.put(Body=str(timestamp))
