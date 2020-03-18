import concurrent.futures
import json
import boto3
import sys
import os
import hashlib


from argparse import ArgumentParser
from shutil import copyfile, rmtree
from urllib.parse import urlparse
from botocore.exceptions import ClientError

sys.path.insert(0, "")  # noqa
from dcp_prototype.backend.wrangling.migrations.utils.util import list_files_in_bucket


s3 = boto3.resource("s3")


def load_files(source_directory, target_directory, checksums={}, clear_if_exists=False):
    """
    Copies and flattens all files (not directories) that exist in `directory` to
    `new_directory_path`. Files will be named based on their uuid and version (with the exception of links.json files
    which are named after their checksum)

    :param source_directory: The directory from which to copy files.
    :param target_directory: The directory to which to copy and de-dupe files.
    :param checksums: A list of checksums of files that have already been processed.
    :param clear_if_exists: If true and the target directory exists, then delete and
    recreate the target directory while removing any existing files (only applicable for non s3 target directories).
    """

    if "s3" in target_directory:
        bucket = urlparse(target_directory).netloc.split(".")[0]
        key = urlparse(target_directory).path[1:]

        print(f"Copying files to S3 bucket and key: {bucket}/{key}")
        checksums = list_files_in_bucket(bucket=bucket, key=key, checksums=True)
        _copy_files_to_s3(source_directory, target_directory, checksums, bucket, key)

    else:
        if not os.path.isdir(target_directory):
            os.mkdir(target_directory)
        elif os.path.isdir(target_directory) and clear_if_exists:
            rmtree(target_directory)

        # Prepopulate checksums with existing file checksums in target directory.
        for filename in os.listdir(target_directory):
            target_file_path = os.path.join(target_directory, filename)
            with open(target_file_path, "rb") as file_reader:
                contents = file_reader.read()
            checksums[hashlib.md5(contents).hexdigest()] = filename

        print(f"Processing directory {source_directory}")
        _copy_files(source_directory, target_directory, checksums)


def _copy_files_to_s3(
    source_directory, target_directory, checksums, bucket, key,
):
    for filename in os.listdir(source_directory):
        # Skip hidden files
        if filename.startswith("."):
            continue

        source_file_path = os.path.join(source_directory, filename)
        if os.path.isdir(source_file_path):
            _copy_files_to_s3(
                source_file_path, target_directory, checksums, bucket, key,
            )
        elif os.path.isfile(source_file_path):
            with open(source_file_path, "rb") as file_reader:
                contents = file_reader.read()
                data = json.loads(contents)
                try:
                    file_checksum = hashlib.md5(contents).hexdigest()
                    uuid = data["provenance"]["document_id"]
                    version = data["provenance"]["update_date"].replace(":", "")
                    file_name = f"{uuid}.{version}.json"
                except KeyError:
                    if source_file_path.split("/")[-1] == "links.json":
                        file_name = file_checksum
                    else:
                        print(f"Issue with file: {source_file_path}")
            if file_checksum not in checksums.values():
                try:
                    s3.meta.client.upload_file(source_file_path, bucket, f"{key}/{file_name}")
                    checksums[file_checksum] = f"{key}.{file_name}"
                except ClientError as e:
                    print(f"error: {e} with file: {file_name}")


def _copy_files(source_directory, target_directory, checksums):
    """
    Recursively copy over files as long as its checksum does not exist in the given list
    of checksums and append the file with the count of the file type as designated by
    the counts in `count_by_file_typ` if the file type exists.
    """

    for filename in os.listdir(source_directory):
        # Skip hidden files
        if filename.startswith("."):
            continue

        source_file_path = os.path.join(source_directory, filename)
        if os.path.isdir(source_file_path):
            print(f"Processing directory {source_file_path}")
            _copy_files(source_file_path, target_directory, checksums)
        elif os.path.isfile(source_file_path):
            with open(source_file_path, "rb") as file_reader:
                contents = file_reader.read()
                data = json.loads(contents)
                uuid = data["provenance"]["document_id"]
                version = data["provenance"]["update_date"].replace(":", "")
                file_name = f"{uuid}.{version}.json"
            file_checksum = hashlib.md5(contents).hexdigest()

            print(f"Processing file: {source_file_path} with checksum {file_checksum}")

            if file_checksum not in checksums:
                checksums["file_checksum"] = f"{uuid}.{version}"
                copyfile(source_file_path, os.path.join(target_directory, file_name))
            else:
                print(f"Skipping over file {filename} because it has already been copied " f"over.")


def load_project_metadata(input_directory, output_directory, project_file, max_workers):
    with open(project_file) as json_file:
        data = json.load(json_file)
        dispatch_executor_class = concurrent.futures.ProcessPoolExecutor
        with dispatch_executor_class(max_workers=max_workers) as executor:
            futures = []
            projects_loaded = 0
            for project in data:
                project = project.replace(" ", "-").replace("/", "-").replace(".", "")
                full_input_directory = f"{input_directory}/{project}/bundles"
                full_output_directory = f"{output_directory}/{project}"
                f = executor.submit(load_files, full_input_directory, full_output_directory)
                futures.append(f)

            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                    projects_loaded += 1
                    print(f"{projects_loaded} projects loaded so far (out of {len(data)})")
                except Exception as e:
                    print(f"Something went wrong: {e}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_directory",
        nargs="+",
        required=True,
        help="A data directory containing the metadata files for one or more projects nested under "
        "{projet_name}/bundles/{bundle_id}/",
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        nargs="+",
        required=True,
        help="An output directory to which the files from the input directory will be "
        "copied to and de-duped. File names will be retained except for the suffix"
        " that is a numbering of the file.",
    )
    parser.add_argument(
        "-t", "--threads", nargs="+", required=False, help="The number of threads for extraction, default is 1",
    )
    parser.add_argument(
        "-p",
        "--project_file",
        nargs="+",
        required=True,
        help="A file containing a list of projects to load metadata for",
    )
    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
    if arguments.output_directory:
        output_directory = arguments.output_directory[0]
    if arguments.threads:
        threads = int(arguments.threads[0])
    else:
        threads = 1
    if arguments.project_file:
        project_file = arguments.project_file[0]

    load_project_metadata(input_directory, output_directory, project_file, threads)
