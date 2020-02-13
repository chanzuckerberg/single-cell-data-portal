from shutil import copyfile, rmtree
import os
import json
import hashlib
from urllib.parse import urlparse
from argparse import ArgumentParser
import boto3
from botocore.exceptions import ClientError

import boto3

s3 = boto3.resource('s3')

COUNT_BY_FILE_TYPE = {
    "analysis_file": 0,
    "analysis_process": 0,
    "analysis_protocol": 0,
    "cell_suspension": 0,
    "dissociation_protocol": 0,
    "donor_organism": 0,
    "enrichment_protocol": 0,
    "library_preparation_protocol": 0,
    "links": 0,
    "process": 0,
    "project": 0,
    "sequence_file": 0,
    "sequencing_protocol": 0,
    "specimen_from_organism": 0,
}


def reset_counts(counts_file_path):
    """
    Resets a dictionary that counts the number of files of that type that the script has
    processed/copied over into a directory that has de-duped them. Every file type will
    have a count of zero.
    """

    global COUNT_BY_FILE_TYPE

    print("Resetting file type counts.")
    for file_type in COUNT_BY_FILE_TYPE.keys():
        COUNT_BY_FILE_TYPE[file_type] = 0
    with open(counts_file_path, "w") as counts_file_pointer:
        counts_file_pointer.write(json.dumps(COUNT_BY_FILE_TYPE))


def delete_directory(directory):
    """
    Delete given directory.
    """

    print(f"Deleting directory: {directory}")
    rmtree(directory)


def squish_files(
    source_directory, target_directory, count_file, checksums=[], clear_if_exists=True, s3=False
):
    """
    Copies and flattens all files (not directories) that exist in `directory` to
    `new_directory_path`. `count_file` is a file containing the last known counts of the
    number of files of each type that has been transferred over in the past to avoid
    overwrites.

    :param source_directory: The directory from which to copy files.
    :param target_directory: The directory to which to copy and de-dupe files.
    :param count_file: A file containing a dictionary that has the count of each file
    processed. This count is used to append a suffix to each file name as it is copied
    over.
    :param checksums: A list of checksums of files that have already been processed.
    :param clear_if_exists: If true and the target directory exists, then delete and
    recreate the target directory while removing any existing files.
    """

    global COUNT_BY_FILE_TYPE

    with open(count_file) as json_file:
        COUNT_BY_FILE_TYPE = json.load(json_file)

    if "s3" in target_directory:
        s3_session = boto3.client("s3")
        bucket = urlparse(target_directory).netloc
        key = urlparse(target_directory).path[1:]
        print(f"Copying files to S3 bucket and key: {bucket}/{key}")
        _copy_files_to_s3(
            source_directory,
            target_directory,
            COUNT_BY_FILE_TYPE,
            checksums,
            s3_session,
            bucket,
            key,
        )

    else:
        if not os.path.isdir(target_directory):
            os.mkdir(target_directory)
        elif os.path.isdir(target_directory) and clear_if_exists:
            rmtree(target_directory)
        elif os.path.isdir(target_directory) and not clear_if_exists:
            for file in os.listdir(target_directory):
                with open(os.path.join(target_directory, file), "rb") as file_reader:
                    contents = file_reader.read()
                checksums.append(hashlib.md5(contents).hexdigest())

        # Prepopulate checksums with existing file checksums in target directory.
        for filename in os.listdir(target_directory):
            target_file_path = os.path.join(target_directory, filename)
            with open(target_file_path, "rb") as file_reader:
                contents = file_reader.read()
            checksums.append(hashlib.md5(contents).hexdigest())

        print(f"Processing directory {source_directory}")
        _copy_files(source_directory, target_directory, COUNT_BY_FILE_TYPE, checksums)

    # Store updated count files
    if count_file:
        json_file_pointer = open(count_file, "w")
        json_file_pointer.write(json.dumps(COUNT_BY_FILE_TYPE))
        json_file_pointer.close()


def _copy_files_to_s3(
    source_directory,
    target_directory,
    count_by_file_type,
    checksums,
    session,
    bucket,
    key,
):
    for filename in os.listdir(source_directory):
        # Skip hidden files
        if filename.startswith("."):
            continue

        source_file_path = os.path.join(source_directory, filename)
        if os.path.isdir(source_file_path):
            print(f"Processing directory {source_file_path}")
            _copy_files_to_s3(
                source_file_path,
                target_directory,
                count_by_file_type,
                checksums,
                session,
                bucket,
                key,
            )
        elif os.path.isfile(source_file_path):
            with open(source_file_path, "rb") as file_reader:
                contents = file_reader.read()
            file_checksum = hashlib.md5(contents).hexdigest()

            print(f"Processing file: {source_file_path} with checksum {file_checksum}")

            if file_checksum not in checksums:
                checksums.append(file_checksum)

                if ".json" in filename:
                    for file_type in COUNT_BY_FILE_TYPE.keys():
                        if file_type in filename:
                            filename = (
                                f"{file_type}_{str(COUNT_BY_FILE_TYPE[file_type])}.json"
                            )
                            COUNT_BY_FILE_TYPE[file_type] += 1

                print(f"Uploading file to S3: {filename}")
                try:
                    # session.upload_file(source_file_path, bucket, key + filename)
                    s3.meta.client.upload_file(source_file_path, 'hca-dcp-one-backup-data', f"{key}/{filename}")

                except ClientError as e:
                    print(e)
            else:
                print(
                    f"Skipping over file {filename} because it has already been copied "
                    f"over."
                )


def _copy_files(source_directory, target_directory, count_by_file_type, checksums):
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
            _copy_files(
                source_file_path, target_directory, count_by_file_type, checksums
            )
        elif os.path.isfile(source_file_path):
            with open(source_file_path, "rb") as file_reader:
                contents = file_reader.read()
            file_checksum = hashlib.md5(contents).hexdigest()

            print(f"Processing file: {source_file_path} with checksum {file_checksum}")

            if file_checksum not in checksums:
                checksums.append(file_checksum)

                if ".json" in filename:
                    for file_type in COUNT_BY_FILE_TYPE.keys():
                        if file_type in filename:
                            filename = (
                                f"{file_type}_{str(COUNT_BY_FILE_TYPE[file_type])}.json"
                            )
                            COUNT_BY_FILE_TYPE[file_type] += 1

                copyfile(source_file_path, os.path.join(target_directory, filename))

            else:
                print(
                    f"Skipping over file {filename} because it has already been copied "
                    f"over."
                )


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
        "-c",
        "--count_file",
        nargs="+",
        required=True,
        help="A count file that keeps track of the number of a particular file that the"
        " script has processed. This file allows for multiple runs of the script "
        "while retaining previous counts.",
    )
    parser.add_argument("-r", "--reset", action="store_true")
    parser.add_argument("-d", "--cleanup", action="store_true")

    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
    if arguments.output_directory:
        output_directory = arguments.output_directory[0]
    if arguments.count_file:
        count_file = arguments.count_file[0]
    should_reset = arguments.reset
    should_delete_source_directory = arguments.cleanup

    if should_reset:
        reset_counts(count_file)
    squish_files(
        input_directory, output_directory, count_file, clear_if_exists=should_reset
    )
    if should_delete_source_directory:
        delete_directory(input_directory)
