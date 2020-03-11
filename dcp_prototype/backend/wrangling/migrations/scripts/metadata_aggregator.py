import os
import queue
import sys
import tarfile
import tempfile
import threading
from urllib.parse import urlparse

import argparse
import boto3
from botocore.exceptions import ClientError

"""This is a script to download all flattened DCP-1.0 metadata files from S3. The files can then be tarred and 
gzipped, then copied back to S3."""


def download_single_file_from_s3(prefix, bucket, filequeue, download_dir):
    """The function runs in a thread.
    While there are files to process, the thread grabs the next filename from the filequeue,
    and downloads that single file into the download_dir.
    """

    session = boto3.session.Session()
    s3_client = session.client("s3")
    while True:
        filename = filequeue.get()

        # check for more files to process
        if filename is None:
            break

        object_name = prefix + "/" + filename
        qsize = filequeue.qsize()
        if qsize % 1000 == 0 and qsize > 0:
            print(f"Working with JSON input file: {object_name}, queued size={qsize}")

        download_dest = os.path.join(download_dir, filename)
        try:
            s3_client.download_file(bucket, object_name, download_dest)
        except Exception as e:
            print(f"Failed {str(e)} : {bucket} {object_name} {download_dest}")

        filequeue.task_done()


def download_all_files_from_s3(s3_uri, num_threads, dirname):
    """Download all the files from s3 for a project into a local directory"""

    bucket_name = urlparse(s3_uri).netloc
    prefix = urlparse(s3_uri).path
    if prefix.startswith("/"):
        prefix = prefix[1:]

    client = boto3.client("s3")
    paginator = client.get_paginator("list_objects")
    page_iterator = paginator.paginate(Bucket=bucket_name, Prefix=prefix)

    filequeue = queue.Queue()
    nfiles = 0
    prefixparts = len(prefix.split("/"))
    for page in page_iterator:
        bucket_objects = page.get("Contents")
        for project_object in bucket_objects:
            object_filename = project_object.get("Key")
            # data_files is a prefix for all the project datafiles (matrix, loom, bam)
            if "data_files" in object_filename:
                continue
            parts = object_filename.split("/")
            # Make sure that only the top level objects are downloaded.
            if len(parts) == prefixparts + 1:
                filequeue.put(parts[-1])
                nfiles += 1
                if nfiles % 1000 == 0:
                    print(nfiles, end=" ", flush=True)

    print(nfiles)
    threads = []
    for _ in range(num_threads):
        thread = threading.Thread(target=download_single_file_from_s3, args=(prefix, bucket_name, filequeue, dirname))
        thread.start()
        threads.append(thread)

    filequeue.join()
    for _ in range(num_threads):
        filequeue.put(None)
    for thread in threads:
        thread.join()

    print("downloaded nfiles=", nfiles)


def create_tarfile(dirname):
    tarfilename = dirname + ".tar.gz"
    nfiles = 0
    with tarfile.open(tarfilename, "w:gz") as tar_handle:
        for filename in os.listdir(dirname):
            tar_handle.add(os.path.join(dirname, filename))
            nfiles += 1

    print("files in tar", nfiles)
    return tarfilename


def upload_to_s3(src_file, dest_s3_uri):
    bucket = urlparse(dest_s3_uri).netloc
    object_name = urlparse(dest_s3_uri).path
    if object_name.startswith("/"):
        object_name = object_name[1:]

    s3_client = boto3.client("s3")
    try:
        print(f"Upload: {src_file}, {bucket}, {object_name}")
        s3_client.upload_file(src_file, bucket, object_name)
        print(f"Upload complete: {object_name}")
    except ClientError as e:
        print("Error:", str(e))
        return False

    return True


def combine_s3_files_into_one_aggregate_file(input_directory, num_threads, dirname, output_object):
    if dirname.endswith("/"):
        dirname = dirname[:-1]
    download_all_files_from_s3(input_directory, num_threads, dirname)
    tarname = create_tarfile(dirname)
    if output_object:
        return upload_to_s3(tarname, output_object)
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input_directory", help="An s3 prefix contains metadata files that were a part of a single DCP 1.0"
    )
    parser.add_argument(
        "-d",
        "--download-dir",
        default=None,
        help="If set, download the meta files into the given directory, else a temp directory will be used",
    )

    parser.add_argument(
        "-o",
        "--output_object",
        default=None,
        help="If set, the tar.gz file of the downloaded data will be upload to S3 at this location",
    )

    parser.add_argument(
        "-t", "--threads", required=False, default=1, type=int, help="Number of threads to use when reading files"
    )

    arguments = parser.parse_args()

    input_directory = arguments.input_directory
    if input_directory.endswith("/"):
        input_directory = input_directory[:-1]

    if arguments.download_dir is None:
        with tempfile.TemporaryDirectory() as dirname:
            ret = combine_s3_files_into_one_aggregate_file(
                input_directory, arguments.threads, dirname, arguments.output_object
            )

    else:
        ret = combine_s3_files_into_one_aggregate_file(
            input_directory, arguments.threads, arguments.download_dir, arguments.output_object
        )

    if ret:
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
