import sys

sys.path.insert(0, "")  # noqa
import boto3
from argparse import ArgumentParser
from urllib.parse import urlparse
import os.path
import json

from dcp_prototype.backend.wrangling.migrations.utils.migration_utils import DatasetMetadata

import queue
import threading
import tarfile


def gather_group_file_list(file_list):
    """
    Process files into groups that can be processed simultaneously. In particular, links need to be processed last.
    This is done so that all entities exist before linkage and linking will not ßoccur between unknown entities.

    All entity files, except link files, look this this:
        UUID.DATE.json
    All link files look like:
        UUID (but without dashes)

    E.g.
    This is an entity file (non-link)
        eeb749e6-ce21-4bf4-a53a-005398a00862.2019-09-20T102839.972Z.json

    This is a link file
        f83094a5b135e7ca4ff9e2974ce36489

    """

    entity_group = [fname for fname in file_list if fname.endswith(".json")]
    link_group = [fname for fname in file_list if not fname.endswith(".json")]
    return [entity_group, link_group]


def consume_file_local(filequeue, dataset_metadata):
    while True:
        filename = filequeue.get()

        # Check for more files to process
        if filename is None:
            break

        qsize = filequeue.qsize()
        if qsize % 1000 == 0 and qsize > 0:
            print(f"Working with JSON input file: {filename}, queued size={qsize}")

        with open(filename) as json_file:
            object_data = json.load(json_file)
            dataset_metadata.add_entity(filename, object_data)
            filequeue.task_done()


def consume_file_tar(filequeue, dataset_metadata):
    while True:
        item = filequeue.get()

        # Check for more files to process
        if item is None:
            break

        filename, buf = item
        qsize = filequeue.qsize()
        if qsize % 1000 == 0 and qsize > 0:
            print(f"Working with tarinfo: queued size={qsize}")

        object_data = json.loads(buf)
        dataset_metadata.add_entity(filename, object_data)
        filequeue.task_done()


def consume_file_s3(prefix, bucket, filequeue, dataset_metadata):
    # Thread local session and client
    session = boto3.session.Session()
    s3_client = session.client("s3")

    while True:
        filename = filequeue.get()

        # Check for more files to process
        if filename is None:
            break

        file_prefix = prefix + filename

        qsize = filequeue.qsize()
        if qsize % 1000 == 0 and qsize > 0:
            print(f"Working with JSON input file: {file_prefix}, queued size={qsize}")
        object = s3_client.get_object(Bucket=bucket, Key=file_prefix)
        object_body = object["Body"].read()
        object_data = json.loads(object_body)

        dataset_metadata.add_entity(filename, object_data)
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

    return dataset_metadata


def generate_metadata_structure_from_dir(input_dir, num_threads, dataset_metadata):
    file_list = os.listdir(input_dir)
    group_file_list = gather_group_file_list(file_list)

    print(f"Files in directory to parse: {len(file_list)}")

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

    return dataset_metadata


def generate_metadata_structure_from_targz(input_source, num_threads, dataset_metadata):
    index = 0
    print("open ", input_source)
    with tarfile.open(input_source, "r:gz") as tar:
        print("get members")
        members = tar.getmembers()
        print("get members done")
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

    for group_files in group_file_list:
        print(f"Files in group to parse: {len(group_files)}")
        index = 0
        for filename, buf in group_files:
            object_data = json.loads(buf)
            dataset_metadata.add_entity(filename, object_data)
            index += 1
            if index % 1000 == 0:
                print("processed ", index)

    return dataset_metadata


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

    parser.add_argument(
        "-o",
        "--output-file",
        default=None,
        help="If provided, the json output will be save at this file.  "
             "If not provided, the json output will go to stdout",
    )

    parser.add_argument(
        "-a",
        "--append",
        action="store_true",
        help="If provided, the json output will be save at this file.  "
             "If not provided, the json output will go to stdout",
    )
    arguments = parser.parse_args()

    input_source = arguments.input_source
    num_threads = arguments.threads

    if "s3" in input_source:

        bucket = urlparse(input_source).netloc
        path = urlparse(input_source).path
        if path.startswith("/"):
            path = path[1:]

        if input_source.endswith(".tar.gz"):
            input_tarfile = os.path.basename(path)
            if os.path.exists(input_tarfile):
                print("Error, cannot overwrite {input_tarfile}")
                sys.exit(1)
            s3_client = boto3.client("s3")
            s3_client.download_file(bucket, path, input_tarfile)
            dataset_metadata = DatasetMetadata()
            generate_metadata_structure_from_targz(input_tarfile, num_threads, dataset_metadata)
        else:
            if not input_source.endswith("/"):
                input_source += "/"

            dataset_metadata = DatasetMetadata()
            generate_metadata_structure_from_s3_uri(input_source, num_threads, dataset_metadata)

    else:
        dataset_metadata = DatasetMetadata()
        if input_source.endswith(".tar.gz"):
            input_tarfile = input_source
            generate_metadata_structure_from_targz(input_tarfile, num_threads, dataset_metadata)
        else:
            generate_metadata_structure_from_dir(input_source, num_threads, dataset_metadata)

    dataset_metadata.process()
    result_project = dataset_metadata.to_dict()

    okay = dataset_metadata.validate(result_project)
    if not okay:
        print("result failed schema validation")
        sys.exit(1)

    for key, value in dataset_metadata.missing.items():
        print("MISSING:", key, value)

    out_dict = result_project
    title = out_dict.get("projects")[0].get("title")
    output_file = arguments.output_file
    if output_file:
        if arguments.append and os.path.exists(output_file):
            with open(output_file) as json_file:
                data = json.load(json_file)
                projects = data.get("projects")
                if projects is None:
                    print(f"expected 'projects' in {output_file}")
                    sys.exit(1)

                do_append = True
                for index, project in enumerate(projects):
                    if project.get("title") == title:
                        print(f"replaced project {title} in {output_file}")
                        projects[index] = out_dict
                        do_append = False
                        break
                if do_append:
                    print(f"append project {title} in {output_file}")
                    projects.append(out_dict)
                out_dict = data
        else:
            print(f"write project {title} to {output_file}")

        with open(arguments.output_file, "w") as json_file:
            json_file.write(json.dumps(out_dict, indent=2))
    else:
        print(json.dumps(out_dict, indent=2))
