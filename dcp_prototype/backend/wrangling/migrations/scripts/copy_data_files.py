import boto3
import concurrent.futures
import json
import os
import typing

from argparse import ArgumentParser

from boto3.s3.transfer import TransferConfig


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


def list_files_in_bucket(project: str, continuation_token: str = None, keys: typing.List[str] = []) -> list:
    """
    Lists all files in a bucket (hard-coded to dunitz-prod-copy) with the prefix {project}/blobs
    :param project: Name of project files belong to
    :param continuation_token: s3 will return a max of 1000 files, pass in the continuation token to get the next 1000 for a project
    :param keys: list of all of the file blob names
    :return: list of all of the file blob names
    """
    s3 = boto3.client("s3")
    if continuation_token:
        response = s3.list_objects_v2(
            Bucket='dunitz-prod-copy',
            Prefix=f"{project}/blobs",
            ContinuationToken=continuation_token
        )
    else:
        response = s3.list_objects_v2(
            Bucket='dunitz-prod-copy',
            Prefix=f"{project}/blobs"
        )
    if response['KeyCount'] > 0:
        for x in response['Contents']:
            keys.append(x['Key'])
    if response['IsTruncated'] is True:
        list_files_in_bucket(project, response['NextContinuationToken'], keys)
    return keys


def copy_between_s3_buckets(key: str, project: str, max_concurrency: int = 20):
    """
    Copy file between s3 buckets
    :param key: File name
    :param project: Project name
    :param max_concurrency: Max concurrent threads for s3 bucket, boto3 default is 10
    """
    config = TransferConfig(max_concurrency=max_concurrency)

    s3 = boto3.resource('s3')
    copy_source = {
        'Bucket': 'org-hca-dss-prod',
        'Key': key
    }
    s3.meta.client.copy(copy_source, 'dunitz-prod-copy', f'{project}/{key}', Config=config)


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
                    project_data_files.append(file_info)
    return project_data_files


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
        for file_info in project_data_files:
            blob_name = compose_blob_key(file_info)
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
        pre_copied_files = list_files_in_bucket(project)

        move_data_files(project_data_files, project, pre_copied_files)
