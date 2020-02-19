import concurrent.futures
import json
import os
import typing

import boto3
from argparse import ArgumentParser


class FileMetadata:
    FORMAT = "format"
    CREATOR_UID = "creator_uid"
    VERSION = "version"
    CONTENT_TYPE = "content-type"
    SIZE = "size"
    CRC32C = "crc32c"
    S3_ETAG = "s3_etag"
    SHA1 = "sha1"
    SHA256 = "sha256"


def compose_blob_key(file_info: typing.Dict[str, str]) -> str:
    """
    Create the key for a blob, given the file metadata.
    :param file_info: This can either be an object that contains the four keys (SHA256, SHA1, S3_ETAG, and CRC32C) in
                      the key_class.
    """
    return "blobs/" + ".".join((
        file_info[FileMetadata.SHA256],
        file_info[FileMetadata.SHA1],
        file_info[FileMetadata.S3_ETAG],
        file_info[FileMetadata.CRC32C]
    ))


def copy_between_s3_buckets(key, project):
    s3 = boto3.resource('s3')
    copy_source = {
        'Bucket': 'org-hca-dss-prod',
        'Key': key
    }
    s3.meta.client.copy(copy_source, 'dunitz-prod-copy', f'{project}/{key}')


def list_data_files_for_project(input_directory, project_name):
    project_data_files = []
    target_directory_path = os.path.join(input_directory, project_name, 'bundle_manifests')
    for file in os.listdir(target_directory_path):
        # skip hidden files
        if file.startswith("."):
            continue
        file_path = os.path.join(target_directory_path, file)
        with open(file_path, "rb") as file_reader:
            contents = file_reader.read()
            files = json.loads(contents)
            for file_info in files['files']:
                if "metadata" in file_info["content-type"]:
                    continue
                else:
                    project_data_files.append(file_info)
    return project_data_files


def move_data_files(project_data_files, project):
    files_moved = 0
    dispatch_executor_class = concurrent.futures.ThreadPoolExecutor
    with dispatch_executor_class(max_workers=1) as executor:
        print(f"{len(project_data_files)} files to copy over for {project}")
        futures = []
        for file_info in project_data_files:
            blob_name = compose_blob_key(file_info)
            f = executor.submit(copy_between_s3_buckets, blob_name, project)
            futures.append(f)

        for future in concurrent.futures.as_completed(futures):
            try:
                extract_result = future.result()
                files_moved += 1
                if files_moved % 50 == 0:
                    print(f"{files_moved} files moved (out of {len(project_data_files)} for project {project}")
            except Exception as e:
                print(f"Something went wrong: {e}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_directory",
        nargs="+",
        required=True,
        help="A data directory containing bundles with files that contain both the "
             "metadata and the data of a single project.",
    )
    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
    for project in os.listdir(input_directory):
        project_data_files = list_data_files_for_project(input_directory, project)
        move_data_files(project_data_files, project)
