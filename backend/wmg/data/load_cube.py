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

# TODO(prathap): The functionality to fully integrate `snapshot_schema_version` will be done
# in a future ticket:
# https://github.com/chanzuckerberg/single-cell-data-portal/issues/5166

###################################### PUBLIC FUNCTIONS #################################


@log_func_runtime
def upload_artifacts_to_s3(
    *, snapshot_source_path: str, snapshot_schema_version: str, snapshot_id: int, is_snapshot_validation_successful=True
) -> str:
    """
    Uploads the generated cubes to S3 and writes metadata that includes
    information about the data schema version, the snapshot id, and the location
    of the snapshot.

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

    sync_command = ["aws", "s3", "sync", snapshot_source_path, snapshot_s3_dest_path]

    subprocess.run(sync_command)

    _write_snapshot_metadata(snapshot_schema_version, snapshot_id, is_snapshot_validation_successful)

    return snapshot_s3_dest_path


###################################### PRIVATE FUNCTIONS #################################


def _get_latest_snapshot_identifier_file_s3_key(snapshot_schema_version: str) -> str:
    """
    Return s3 key of latest_snapshot_identifier file
    """
    data_root_path = _get_wmg_s3_root_dir_path()
    _get_wmg_s3_data_schema_dir_path(snapshot_schema_version)
    file_name = "latest_snapshot_identifier"

    # TODO(prathap) - change to: return f"{data_root_path}/{data_schema_dir_path}/{file_name}"
    # when working on: https://github.com/chanzuckerberg/single-cell-data-portal/issues/5166
    return f"{data_root_path}/{file_name}"


def _get_latest_snapshot_run_file_s3_key() -> str:
    """
    Return s3 key of latest_snapshot_run file
    """
    data_root_path = _get_wmg_s3_root_dir_path()
    file_name = "latest_snapshot_run"

    return f"{data_root_path}/{file_name}"


def _get_wmg_s3_bucket_uri():
    """
    Get s3 bucket uri
    """
    return f"s3://{wmg_bucket_name}"


def _get_wmg_s3_root_dir_path():
    """
    Get root directory in the s3 bucket under which the data artifact directories will be written to
    """
    return stack_name.strip("/") if stack_name else ""


def _get_wmg_s3_data_schema_dir_path(snapshot_schema_version: str) -> str:
    """
    Get path for a particular data schema version under which data should be located
    """
    return f"snapshots/{snapshot_schema_version}"


def _get_wmg_snapshot_s3_path(
    snapshot_schema_version: str, snapshot_id: str, is_snapshot_validation_successful: bool
) -> str:
    """
    Get the path of the directory in which the data artifact is located
    """
    snapshot_dir = snapshot_id if is_snapshot_validation_successful else f"{snapshot_id}_validation_failed_snapshot"

    data_root_path = _get_wmg_s3_root_dir_path()
    _get_wmg_s3_data_schema_dir_path(snapshot_schema_version)

    # TODO(prathap) - change to: return f"{data_root_path}/{data_schema_dir_path}/{snapshot_dir}"
    # when working on: https://github.com/chanzuckerberg/single-cell-data-portal/issues/5166
    return f"{data_root_path}/{snapshot_dir}"


def _get_wmg_snapshot_s3_fullpath(
    snapshot_schema_version: str, snapshot_id: str, is_snapshot_validation_successful: bool
) -> str:
    """
    Get the fully qualified s3 location of the data artifact
    """
    snapshot_path = _get_wmg_snapshot_s3_path(snapshot_schema_version, snapshot_id, is_snapshot_validation_successful)

    return f"{_get_wmg_s3_bucket_uri()}/{snapshot_path}"


def _make_snapshot_active(snapshot_schema_version, snapshot_id):
    """
    Update the latest snapshot reference file to point to the newly generated snapshot and remove older snapshots
    """
    _write_latest_snapshot_identifier_file(snapshot_schema_version, snapshot_id)

    data_root_path = _get_wmg_s3_root_dir_path()
    _get_wmg_s3_data_schema_dir_path(snapshot_schema_version)

    # TODO(prathap): Change to - s3_key_prefix = f"{data_root_path}/{data_schema_dir_path}"
    # when working on: https://github.com/chanzuckerberg/single-cell-data-portal/issues/5166
    s3_key_prefix = data_root_path

    _remove_oldest_datasets(s3_key_prefix)


def _remove_oldest_datasets(s3_key_prefix):
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
    key_path = _get_latest_snapshot_identifier_file_s3_key(snapshot_schema_version)

    _write_value_to_s3_key(key_path=key_path, value=snapshot_id)


def _write_latest_snapshot_run_file(
    snapshot_schema_version: str, snapshot_id: str, is_snapshot_validation_successful: str
):
    """
    Update latest_snapshot_run file with the path to the latest snapshot folder generated
    """
    key_path = _get_latest_snapshot_run_file_s3_key()
    snapshot_path = _get_wmg_snapshot_s3_path(snapshot_schema_version, snapshot_id, is_snapshot_validation_successful)

    _write_value_to_s3_key(key_path=key_path, value=snapshot_path)


def _write_snapshot_metadata(snapshot_schema_version: str, snapshot_id: str, is_snapshot_validation_successful: str):
    """
    Write metadata about the data artifact generated
    """
    if is_snapshot_validation_successful:
        _make_snapshot_active(snapshot_schema_version, snapshot_id)

    _write_latest_snapshot_run_file(snapshot_schema_version, snapshot_id, is_snapshot_validation_successful)
