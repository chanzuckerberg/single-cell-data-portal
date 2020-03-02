import sys

sys.path.insert(0, "")  # noqa
import boto3
from argparse import ArgumentParser
from pandas.io.json import json_normalize
from urllib.parse import urlparse
import os.path
from flatten_json import flatten
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_dataset_metadata import (  # noqa
    OldDatasetMetadata,
)
import json
from dcp_prototype.backend.ledger.code.common.ledger_orm import DBSessionMaker
from dcp_prototype.backend.wrangling.migrations.utils.constants import ENTITY_TYPES

import queue
import threading
import time
import tarfile


def get_entity_type(object_name, object_data):
    desc = object_data.get("describedBy")
    if desc:
        path = urlparse(desc).path
        entity_type = os.path.basename(path)
        if entity_type not in ENTITY_TYPES:
            if entity_type in (
                "enrichment_protocol",
                "analysis_process",
                "analysis_protocol",
                "process",
                "dissociation_protocol",
            ):
                return None
            raise RuntimeError(f"object {object_name} did not have an expected describedBy field: {entity_type}")
        return entity_type
    else:
        raise RuntimeError(f"object {object_name} did not have a describedBy field")


def gather_group_file_list(file_list):
    """
    process files into groups that can be processed simultaneously.
    In particular, links need to be processed last.
    This is done so that all entities exist before linkage and linking will not
    occur between unknown entities.

    All entity files, except link files, look this this:
        UUID.DATE.json
    All link files look like:
        UUID (but without dashes)

    E.g.
    this is an entity file (non-link)
        eeb749e6-ce21-4bf4-a53a-005398a00862.2019-09-20T102839.972Z.json

    this is a link file
        f83094a5b135e7ca4ff9e2974ce36489

    """

    entity_group = [fname for fname in file_list if fname.endswith(".json")]
    link_group = [fname for fname in file_list if not fname.endswith(".json")]
    return [entity_group, link_group]


def process_file_data(filename, object_data, dataset_metadata):
    tsv_generated_data_frame = json_normalize(flatten(object_data, separator="."))
    entity_type = get_entity_type(filename, object_data)
    if entity_type:
        for row in tsv_generated_data_frame.iterrows():
            metadata_row = row[1]
            dataset_metadata.parse_flattened_row_of_json(metadata_row, entity_type)


def consume_file_local(filequeue, dataset_metadata):
    while True:
        filename = filequeue.get()

        # check for more files to process
        if filename is None:
            break

        qsize = filequeue.qsize()
        if qsize % 1000 == 0 and qsize > 0:
            print(f"Working with JSON input file: {filename}, queued size={qsize}")

        with open(filename) as json_file:
            object_data = json.load(json_file)
            process_file_data(filename, object_data, dataset_metadata)
            filequeue.task_done()


def consume_file_tar(filequeue, dataset_metadata):
    while True:
        item = filequeue.get()

        # check for more files to process
        if item is None:
            break

        filename, buf = item
        qsize = filequeue.qsize()
        if qsize % 1000 == 0 and qsize > 0:
            print(f"Working with tarinfo: queued size={qsize}")

        object_data = json.loads(buf)
        process_file_data(filename, object_data, dataset_metadata)
        filequeue.task_done()


def consume_file_s3(prefix, bucket, filequeue, dataset_metadata):

    # thread local session and client
    session = boto3.session.Session()
    s3_client = session.client("s3")

    while True:
        filename = filequeue.get()

        # check for more files to process
        if filename is None:
            break

        file_prefix = prefix + filename

        qsize = filequeue.qsize()
        if qsize % 1000 == 0 and qsize > 0:
            print(f"Working with JSON input file: {file_prefix}, queued size={qsize}")
        object = s3_client.get_object(Bucket=bucket, Key=file_prefix)
        object_body = object["Body"].read()
        object_data = json.loads(object_body)

        process_file_data(filename, object_data, dataset_metadata)
        filequeue.task_done()


def generate_metadata_structure_from_s3_uri(s3_uri, num_threads, dataset_metadata):
    """
    Transforms a project's metadata into the DCP 2.0 metadata schema. Metadata
    exists primarily in JSON files in the given S3 bucket.
    """

    bucket_name = urlparse(s3_uri).netloc
    PREFIX = urlparse(s3_uri).path
    if PREFIX.startswith("/"):
        PREFIX = PREFIX[1:]

    s3_client = boto3.client("s3")
    paginator = s3_client.get_paginator("list_objects")
    page_iterator = paginator.paginate(Bucket=bucket_name, Prefix=PREFIX)

    file_list = []
    print("Gather objects: ", end="", flush=True)
    for page in page_iterator:
        bucket_objects = page.get("Contents")
        for object in bucket_objects:
            object_filename = object.get("Key")
            file_list.append(object_filename.split("/")[-1])
            if len(file_list) % 1000 == 0:
                print(len(file_list), end=" ", flush=True)

    print()
    group_file_list = gather_group_file_list(file_list)

    print(f"Files in directory to parse: {len(file_list)}")

    tstart = time.time()
    for group_files in group_file_list:
        print(f"Files in group to parse: {len(group_files)}")
        filequeue = queue.Queue()

        for filename in group_files:
            filequeue.put(filename)

        threads = []
        for _ in range(num_threads):
            thread = threading.Thread(target=consume_file_s3, args=(PREFIX, bucket_name, filequeue, dataset_metadata))
            thread.start()
            threads.append(thread)

        filequeue.join()
        for _ in range(num_threads):
            filequeue.put(None)
        for thread in threads:
            thread.join()

    tend = time.time()
    print(f"Process file time (t={num_threads}): {(tend-tstart)}")
    return dataset_metadata


def generate_metadata_structure_from_dir(input_dir, num_threads, dataset_metadata):

    file_list = os.listdir(input_dir)
    group_file_list = gather_group_file_list(file_list)

    print(f"Files in directory to parse: {len(file_list)}")

    tstart = time.time()
    for group_files in group_file_list:
        print(f"Files in group to parse: {len(group_files)}")
        filequeue = queue.Queue()

        for filename in group_files:
            filequeue.put(os.path.join(input_dir, filename))

        threads = []
        for _ in range(num_threads):
            thread = threading.Thread(target=consume_file_local, args=(filequeue, dataset_metadata))
            thread.start()
            threads.append(thread)

        filequeue.join()
        for _ in range(num_threads):
            filequeue.put(None)
        for thread in threads:
            thread.join()

    tend = time.time()
    print(f"Process file time (t={num_threads}): {(tend - tstart)}")
    return dataset_metadata


def generate_metadata_structure_from_targz(input_source, num_threads, dataset_metadata):

    index = 0
    with tarfile.open(input_source, "r:gz") as tar:
        members = tar.getmembers()
        entity_group = []
        link_group = []
        for index, member in enumerate(members):
            buf = tar.extractfile(member).read()
            if member.name.endswith(".json"):
                entity_group.append((member.name, buf))
            else:
                link_group.append((member.name, buf))

            if index > 0 and index % 1000 == 0:
                print(index, end=" ", flush=True)

    print(index)
    group_file_list = [entity_group, link_group]
    print(f"Files in directory to parse: {len(members)}")

    tstart = time.time()
    for group_files in group_file_list:
        print(f"Files in group to parse: {len(group_files)}")
        index = 0
        for filename, buf in group_files:
            object_data = json.loads(buf)
            process_file_data(filename, object_data, dataset_metadata)
            index += 1
            if index % 1000 == 0:
                print("processed ", index)

    tend = time.time()
    print(f"Process file time (t={num_threads}): {(tend - tstart)}")
    return dataset_metadata


def export_old_metadata_to_s3_orm(dataset_metadata: OldDatasetMetadata):
    # Upload to DB
    print(f"Committing dataset to Ledger DB")
    session_maker = DBSessionMaker()
    dataset_metadata.export_to_database(session_maker)
    print(f"Completed updating Ledger DB with dataset metadata")


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input-source",
        required=True,
        help="input_source contains all metadata files that were part of a single DCP 1.0 project."
        " It may be one of the following:  An s3 bucket containing json files."
        " A directory containing json files."
        " An s3 object in the tar.gz format containing the json files."
        " Or a local tar.gz file containing the json files.",
    )

    parser.add_argument(
        "-t", "--threads", required=False, default=1, type=int, help="Number of threads to use when reading files"
    )

    parser.add_argument("-s", "--s3uri", help="The source of the data (s3uri)")

    arguments = parser.parse_args()

    input_source = arguments.input_source
    input_tarfile = None
    num_threads = arguments.threads

    s3_uri = arguments.s3uri

    if "s3" in input_source:

        bucket = urlparse(input_source).netloc
        path = urlparse(input_source).path
        if path.startswith("/"):
            path = path[1:]

        if input_source.endswith(".tar.gz"):
            if s3_uri is None:
                dpath = path.replace(".tar.gz", "")
                s3_uri = f"s3://{bucket}/{dpath}"

            input_tarfile = os.path.basename(path)
            if os.path.exists(input_tarfile):
                print("Error, cannot overwrite {input_tarfile}")
                sys.exit(1)
            s3_client = boto3.client("s3")
            s3_client.download_file(bucket, path, input_tarfile)
            if s3_uri is None:
                print("You must specify an s3uri")
                sys.exit(1)
            dataset_metadata = OldDatasetMetadata(s3_uri=s3_uri)
            generate_metadata_structure_from_targz(input_tarfile, num_threads, dataset_metadata)
        else:
            if not input_source.endswith("/"):
                input_source += "/"

            if s3_uri is None:
                s3_uri = f"s3://{bucket}/{path}"

            dataset_metadata = OldDatasetMetadata(s3_uri=s3_uri)
            generate_metadata_structure_from_s3_uri(input_source, num_threads, dataset_metadata)

    else:
        if s3_uri is None:
            print("You must specify an s3uri")
            sys.exit(1)
        dataset_metadata = OldDatasetMetadata(s3_uri=s3_uri)
        if input_source.endswith(".tar.gz"):
            input_tarfile = input_source
            generate_metadata_structure_from_targz(input_tarfile, num_threads, dataset_metadata)
        else:
            generate_metadata_structure_from_dir(input_source, num_threads, dataset_metadata)

    tstart = time.time()
    export_old_metadata_to_s3_orm(dataset_metadata)
    tend = time.time()
    print("export_old_metadata_to_s3_orm:", (tend - tstart))


if __name__ == "__main__":
    main()
