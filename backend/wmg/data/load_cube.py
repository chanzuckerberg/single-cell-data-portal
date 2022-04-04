import boto3
import os
import subprocess


stack_name = os.environ.get("REMOTE_DEV_PREFIX")
wmg_bucket_name = os.environ.get("WMG_BUCKET")
artifact_bucket_name = os.environ.get("ARTIFACT_BUCKET")


def _get_wmg_bucket_path():
    if stack_name:
        return f"s3://{wmg_bucket_name}{stack_name}"
    else:
        return f"s3://{wmg_bucket_name}"


def update_s3_resources(snapshot_path, timestamp):
    """
    Copy cube and cube related files to s3 under the current timestamp
    """
    upload_artifacts_to_s3(snapshot_path, timestamp)
    update_latest_snapshot_identifier(timestamp)
    remove_oldest_datasets(timestamp)


def remove_oldest_datasets(timestamp):
    """
    Remove all data stored under the oldest timestamp
    """
    s3 = boto3.resource("s3")
    wmg_bucket = s3.Bucket(wmg_bucket_name)
    objects = wmg_bucket.objects.filter(Prefix=stack_name.strip("/")) if stack_name else wmg_bucket.objects.all()

    def enumerate_timestamps_and_objects(objects):
        for obj in objects:
            try:
                tokens = obj.key.split("/")
                timestamp_token = tokens[1] if stack_name else tokens[0]
                is_timestamp = timestamp_token[:10].isdigit()
                if is_timestamp:
                    yield (timestamp_token, obj)
            except Exception:
                pass

    candidate_to_delete = [obj for obj in enumerate_timestamps_and_objects(objects)]

    timestamps = sorted(list(set([x[0] for x in candidate_to_delete])))

    if len(timestamps) > 2:
        timestamps_to_delete = list(timestamps)[:-2]
    else:
        timestamps_to_delete = []

    for timestamp, object in candidate_to_delete:
        if timestamp in timestamps_to_delete:
            object.delete()


def upload_artifacts_to_s3(snapshot_path, timestamp):
    sync_command = ["aws", "s3", "sync", snapshot_path, f"{_get_wmg_bucket_path()}/{timestamp}"]
    subprocess.run(sync_command)


def update_latest_snapshot_identifier(timestamp):
    """
    Update static timestamp in s3 to match the latest
    """
    s3 = boto3.resource("s3")
    file_key = "latest_snapshot_identifier"
    file_path = f"{stack_name.strip('/')}/{file_key}" if stack_name else f"{file_key}"
    object = s3.Object(wmg_bucket_name, file_path)
    object.put(Body=str(timestamp))
