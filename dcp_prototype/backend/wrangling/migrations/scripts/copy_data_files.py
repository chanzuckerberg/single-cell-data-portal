import sys

import boto3
import concurrent.futures
import json
import os
import typing

from argparse import ArgumentParser

from boto3.s3.transfer import TransferConfig
sys.path.insert(0, "")  # noqa
from dcp_prototype.backend.wrangling.migrations.utils.util import list_files_in_bucket

s3 = boto3.resource('s3')
SOURCE_BUCKET = 'org-hca-dss-prod'
DESTINATION_BUCKET = 'dunitz-prod-copy'


class FileMetadata:
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


def copy_between_s3_buckets(key: str, project: str, max_concurrency: int = 20):
    """
    Copy file between s3 buckets
    :param key: File name
    :param project: Project name
    :param max_concurrency: Max concurrent threads for s3 bucket, boto3 default is 10
    """
    config = TransferConfig(max_concurrency=max_concurrency)

    copy_source = {
        'Bucket': SOURCE_BUCKET,
        'Key': key
    }
    s3.meta.client.copy(copy_source, DESTINATION_BUCKET, f'{project}/{key}', Config=config)


def list_data_files_for_project(directory_path: str) -> typing.List[typing.Dict[str, str]]:
    """
    Read all bundle manifests in a particular directory and return a list of non metadata files, file information
    :param directory_path: Path to bundle manifest files
    :return: List of dictionaries containing file information (hash, size, ids etc.) for all non metadata files for that project
    """
    project_data_files = []
    for file in os.listdir(directory_path):
        # skip hidden files
        if file.startswith("."):
            continue
        file_path = os.path.join(directory_path, file)
        with open(file_path, "rb") as file_reader:
            contents = file_reader.read()
            files = json.loads(contents)
            for file_info in files['files']:
                if "metadata" in file_info["content-type"]:
                    continue
                else:
                    project_data_files.append(compose_blob_key(file_info))
    return list(set(project_data_files))


def move_data_files(project_data_files: typing.List[str], project: str, precopied_files: typing.List[str]):
    """
    Parallelize calls to copy files between s3 buckets (if not already in destination bucket)
    :param project_data_files: List of files to copy
    :param project: Name of project (for destination key)
    :param precopied_files: List of files already in destination bucket
    """
    files_moved = 0
    already_there = 0
    dispatch_executor_class = concurrent.futures.ThreadPoolExecutor
    with dispatch_executor_class(max_workers=20) as executor:
        print(f"{len(project_data_files)} files to copy over for {project}")
        futures = []
        for blob_name in project_data_files:
            key = f"{project}/{blob_name}"
            if key in precopied_files:
                already_there += 1
                continue
            f = executor.submit(copy_between_s3_buckets, blob_name, project)
            futures.append(f)

        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
                files_moved += 1
                if files_moved % 500 == 0:
                    print(
                        f"{files_moved} (of {len(project_data_files)}) files moved for project {project}. *{already_there} files were already in the bucket")
            except Exception as e:
                print(f"Something went wrong: {e}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_directory",
        nargs="+",
        required=True,
        help="A data directory containing multiple projects with bundle manifest file(s) "
             "nested under {project_name}/bundle_manifests/",
    )
    arguments = parser.parse_args()
    if arguments.input_directory:
        input_directory = arguments.input_directory[0]

    for project in os.listdir(input_directory):
        directory_path = os.path.join(input_directory, project, 'bundle_manifests')
        project_data_files = list_data_files_for_project(directory_path)
        pre_copied_files = list_files_in_bucket(DESTINATION_BUCKET, f"{project}/blobs")
        move_data_files(project_data_files, project, pre_copied_files)
