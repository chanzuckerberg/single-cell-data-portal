import logging
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

# root directory in the s3 bucket under which the data artifact directories will be written to
wmg_s3_root_dir_path = stack_name.strip("/") if stack_name else ""

wmg_bucket_name = os.environ.get("WMG_BUCKET")
wmg_s3_bucket_uri = f"s3://{wmg_bucket_name}"

logger = logging.getLogger(__name__)

###################################### PUBLIC FUNCTIONS #################################


@log_func_runtime
def upload_artifacts_to_s3(
    *, snapshot_source_path: str, snapshot_schema_version: str, snapshot_id: int, is_snapshot_validation_successful=True
) -> str:
    """
    Uploads the generated cubes to S3 and writes metadata that includes
    information about the data schema version, the snapshot id, and the location
    of the snapshot.

    The <destination_path> is set as s3://<environment_bucket_name>/snapshots/<data_schema_version>/

    1. Snapshot artifacts are stored in <destination_path>/<snapshot_id>/
       A. Note that, as of now, <snapshot_id> is a integer representing the time elapsed
          in seconds since epoch
       B. If the snapshot generation pipeline is run in a "remote dev environment", then the snapshot
          is stored in s3://<environment_bucket_name>/<remote_stack_name>/snapshots/<data_schema_version>/<snapshot_id>/
       C. If validation checks on the snapshot fails, the data artifact will be stored in
          s3://<environment_bucket_name>/snapshots/<data_schema_version>/<snapshot_id>_validation_failed_snapshot/
          (or s3://<environment_bucket_name>/<remote_stack_name>/<snapshot_id>_validation_failed_snapshot/)

    2. The file, s3://<environment_bucket_name>/snapshots/<data_schema_version>/latest_snapshot_identifer, contains
       the latest <snapshot_id> that has passed validation for that schema version. If snapshot generation pipeline
       is run in a "remote dev environment" then the file will be at path:
       s3://<environment_bucket_name>/<remote_stack_name>/snapshots/<data_schema_version>/latest_snapshot_identifer

       NOTE that `latest_snapshot_identifier` will not be updated if validation checks fails on a newly generated
       data artifact

    3. The file, s3://<environment_bucket_name>/latest_snapshot_run, contains the path to the folder containing the
        latest data generated (whether the generated data passed validation or not)

    4. After each generation of the data artifact and successful validation checks, all but the two latest
       snapshots are deleted

    Parameters
    ----------
    snapshot_source_path: The current source path of the cubes that need to be uploaded to s3
    snapshot_schema_version: The schema version of the cube data
    snapshot_id: The snapshot id is a timestamp in seconds
    is_snapshot_validation_successful: Boolean indicating whether data validation ran sucessfully

    Returns
    -------
    snapshot_s3_dest_path: A string that contains the location of the uploaded data
    """
    snapshot_id = str(snapshot_id)

    # cleanup files we do not wish to upload
    shutil.rmtree(f"{snapshot_source_path}/{INTEGRATED_ARRAY_NAME}")
    shutil.rmtree(f"{snapshot_source_path}/{OBS_ARRAY_NAME}")
    shutil.rmtree(f"{snapshot_source_path}/{VAR_ARRAY_NAME}")

    snapshot_s3_dest_path = _get_wmg_snapshot_s3_fullpath(
        snapshot_schema_version, snapshot_id, is_snapshot_validation_successful
    )

    logger.info(f"Writing snapshot data to {snapshot_s3_dest_path}")
    sync_command = ["aws", "s3", "sync", snapshot_source_path, snapshot_s3_dest_path, "--quiet"]
    subprocess.run(sync_command)

    _write_snapshot_metadata(snapshot_schema_version, snapshot_id, is_snapshot_validation_successful)

    return snapshot_s3_dest_path


###################################### PRIVATE FUNCTIONS #################################


def _get_wmg_s3_data_schema_dir_path(snapshot_schema_version: str) -> str:
    """
    Get path for a particular data schema version under which data should be located
    """
    data_schema_dir_path = f"snapshots/{snapshot_schema_version}"

    if wmg_s3_root_dir_path:
        data_schema_dir_path = f"{wmg_s3_root_dir_path}/{data_schema_dir_path}"

    return data_schema_dir_path


def _get_wmg_snapshot_s3_path(
    snapshot_schema_version: str, snapshot_id: str, is_snapshot_validation_successful: bool
) -> str:
    """
    Get the path of the directory in which the data artifact is located
    """
    snapshot_dir = snapshot_id if is_snapshot_validation_successful else f"{snapshot_id}_validation_failed_snapshot"

    data_schema_dir_path = _get_wmg_s3_data_schema_dir_path(snapshot_schema_version)

    return f"{data_schema_dir_path}/{snapshot_dir}" if data_schema_dir_path else snapshot_dir


def _get_wmg_snapshot_s3_fullpath(
    snapshot_schema_version: str, snapshot_id: str, is_snapshot_validation_successful: bool
) -> str:
    """
    Get the fully qualified s3 location of the data artifact
    """
    snapshot_path = _get_wmg_snapshot_s3_path(snapshot_schema_version, snapshot_id, is_snapshot_validation_successful)

    return f"{wmg_s3_bucket_uri}/{snapshot_path}"


def _make_snapshot_active(snapshot_schema_version, snapshot_id):
    """
    Update the latest snapshot reference file to point to the newly generated snapshot and remove older snapshots
    """
    _write_latest_snapshot_identifier_file(snapshot_schema_version, snapshot_id)

    data_schema_dir_path = _get_wmg_s3_data_schema_dir_path(snapshot_schema_version)

    _remove_oldest_datasets(s3_key_prefix=data_schema_dir_path)


def _remove_oldest_datasets(*, s3_key_prefix):
    """
    Remove all snapshots that are older the 2 most recent snapshots
    """
    s3 = boto3.resource("s3")
    wmg_bucket = s3.Bucket(wmg_bucket_name)
    objects = wmg_bucket.objects.filter(Prefix=s3_key_prefix.strip("/")) if s3_key_prefix else wmg_bucket.objects.all()

    # If there is a s3_key_prefix = "x/y/z" and key = "x/y/z/k"
    # then splitting by "/" gives ['x', 'y', 'z', 'k'] which
    # isolates "k" in index 3 = len("x/y/z".split('/'))
    idx = len(s3_key_prefix.split("/")) if s3_key_prefix else 0

    candidate_to_delete = []
    for obj in objects:
        try:
            tokens = obj.key.split("/")
            timestamp_token = tokens[idx]
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


def _write_value_to_s3_key(*, key_path: str, value: str):
    """
    Write a value to an s3 key
    """
    s3 = boto3.resource("s3")
    s3_object = s3.Object(wmg_bucket_name, key_path)
    s3_object.put(Body=value)


def _write_latest_snapshot_identifier_file(snapshot_schema_version: str, snapshot_id: str):
    """
    Update latest_snapshot_identifier file with latest snapshot_id
    """
    data_schema_dir_path = _get_wmg_s3_data_schema_dir_path(snapshot_schema_version)
    file_name = "latest_snapshot_identifier"

    key_path = f"{data_schema_dir_path}/{file_name}" if data_schema_dir_path else file_name

    logger.info(f"Writing latest snapshot identifier file to {key_path} with value {snapshot_id}")
    _write_value_to_s3_key(key_path=key_path, value=snapshot_id)


def _write_latest_snapshot_run_file(
    snapshot_schema_version: str, snapshot_id: str, is_snapshot_validation_successful: str
):
    """
    Update latest_snapshot_run file with the path to the latest snapshot folder generated
    """
    file_name = "latest_snapshot_run"

    key_path = f"{wmg_s3_root_dir_path}/{file_name}" if wmg_s3_root_dir_path else file_name
    snapshot_path = _get_wmg_snapshot_s3_path(snapshot_schema_version, snapshot_id, is_snapshot_validation_successful)

    logger.info(f"Writing latest snapshot run file to {key_path} with value {snapshot_path}")
    _write_value_to_s3_key(key_path=key_path, value=snapshot_path)


def _write_snapshot_metadata(snapshot_schema_version: str, snapshot_id: str, is_snapshot_validation_successful: str):
    """
    Write metadata about the data artifact generated
    """
    if is_snapshot_validation_successful:
        logger.info(f"Snapshot {snapshot_id} passed validation, marking as active")
        _make_snapshot_active(snapshot_schema_version, snapshot_id)
    else:
        logger.info(f"Snapshot {snapshot_id} failed validation")

    _write_latest_snapshot_run_file(snapshot_schema_version, snapshot_id, is_snapshot_validation_successful)
