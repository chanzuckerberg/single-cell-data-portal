import sys

sys.path.insert(0, "")  # noqa
import boto3
from argparse import ArgumentParser
from pandas.io.json import json_normalize
from os import listdir
from urllib.parse import urlparse
import os.path
from botocore.exceptions import ClientError
from flatten_json import flatten
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_dataset_metadata import (  # noqa
    OldDatasetMetadata,
)
import json
from dcp_prototype.backend.ledger.code.common.ledger_orm import DBSessionMaker

BUCKET_NAME = "hca-dcp-one-backup-data"
PREFIX = "single_cell_transcriptome_analysis_of_human_pancreas_metadata"

ENTITY_TYPES = [
    "donor_organism",
    "specimen_from_organism",
    "cell_suspension",
    "library_preparation_protocol",
    "project",
    "sequence_file",
    "links",
    "sequencing_protocol",
    "analysis_file",
]

S3_CLIENT = boto3.client("s3")
S3_RESOURCE = boto3.resource("s3")


def generate_metadata_tsv_name(project_name):
    """
    Given a project name, generates a name for the spreadsheet that will be generated
    based on this project's metadata.
    """
    replace_spaces = project_name.replace(" ", "_")
    replace_slashes = replace_spaces.replace("/", "_")
    return replace_slashes + "_metadata.xlsx"


def get_entity_type(filename):
    for key in ENTITY_TYPES:
        if key in filename:
            return key


def order_file_list(file_list):
    """
    Order files to process to construct metadata schema such that links are processed
    last. This is done so that all entities exist before linkage and linking will not
    occur between unknown entities.
    """
    ordered_file_list = []

    for file in file_list:
        if "donor_organism" in file:
            ordered_file_list.append(file)
    for file in file_list:
        if "specimen" in file:
            ordered_file_list.append(file)
    for file in file_list:
        if "cell_suspension" in file:
            ordered_file_list.append(file)
    for file in file_list:
        if file not in ordered_file_list and "links" not in file:
            ordered_file_list.append(file)
    for file in file_list:
        if file not in ordered_file_list and "links" in file:
            ordered_file_list.append(file)

    return ordered_file_list


def generate_metadata_structure_from_s3_uri(s3_uri):
    """
    Transforms a project's metadata into the DCP 2.0 metadata schema. Metadata
    exists primarily in JSON files in the given S3 bucket.
    """

    dataset_metadata = OldDatasetMetadata(s3_uri=s3_uri)

    BUCKET_NAME = urlparse(s3_uri).netloc
    PREFIX = urlparse(s3_uri).path
    if PREFIX.startswith("/"):
        PREFIX = PREFIX[1:]

    paginator = S3_CLIENT.get_paginator("list_objects")
    page_iterator = paginator.paginate(Bucket=BUCKET_NAME, Prefix=PREFIX)

    file_list = []
    for page in page_iterator:
        bucket_objects = page.get("Contents")
        print(len(bucket_objects))

        for object in bucket_objects:
            object_filename = object.get("Key")
            if should_process_file(object_filename):
                file_list.append(object_filename.split("/")[-1])

    file_list = order_file_list(file_list)

    print(f"Files in directory to parse: {file_list}")

    for filename in file_list:
        if not should_process_file(filename):
            continue

        file_prefix = PREFIX + filename

        print(f"Working with JSON input file: {file_prefix}")
        object = S3_RESOURCE.Object(BUCKET_NAME, file_prefix)
        object_body = object.get()["Body"].read()
        tsv_generated_data_frame = json_normalize(
            flatten(json.loads(object_body), separator=".")
        )
        for row in tsv_generated_data_frame.iterrows():
            metadata_row = row[1]
            dataset_metadata.parse_flattened_row_of_json(
                metadata_row, get_entity_type(filename)
            )

    return dataset_metadata


def generate_metadata_structure(input_directory):
    """
    Transforms a project's metadata into the DCP 2.0 metadata schema. Metadata exists
    primarily in JSON files in the given input directory.
    """

    if "s3" in input_directory:
        return generate_metadata_structure_from_s3_uri(input_directory)

    dataset_metadata = OldDatasetMetadata(s3_uri=f"s3://{BUCKET_NAME}/{PREFIX}")

    ordered_files = order_file_list(listdir(input_directory))
    print(f"Files in directory to parse: {ordered_files}")

    for file in ordered_files:
        full_file_path = os.path.join(input_directory, file)
        if not should_process_file(full_file_path):
            continue
        if os.path.isfile(full_file_path):
            print(f"Working with JSON input file at {full_file_path}")

            with open(full_file_path) as json_file_object:
                input_file_contents = json_file_object.read()
            tsv_generated_data_frame = json_normalize(
                flatten(json.loads(input_file_contents), separator=".")
            )

            for row in tsv_generated_data_frame.iterrows():
                metadata_row = row[1]
                dataset_metadata.parse_flattened_row_of_json(
                    metadata_row, get_entity_type(full_file_path)
                )

    dataset_metadata.save()

    # print(f"Processing dataset: {dataset_metadata.project_name}")
    # dataset_metadata.export_to_spreadsheet(
    #    generate_metadata_tsv_name(dataset_metadata.project_name)
    # )
    print(
        f"Completed generating spreadsheet: "
        f"{generate_metadata_tsv_name(dataset_metadata.project_name)}"
    )

    return dataset_metadata


def should_process_file(full_file_path):
    """
    Returns whether a file, originally from a bundle, needs to be processed in order to
    construct the metadata schema representing the project.
    """
    if "analysis_p" in full_file_path:
        return False
    if "zarr" in full_file_path:
        return False
    if ".json" not in full_file_path:
        return False
    return True


def export_old_metadata_to_s3_orm(
    dataset_metadata: OldDatasetMetadata, input_directory, should_upload_to_s3
):
    """
    Uploads the dataset metadata schema that has been constructed directly to the Ledger
    database using its ORM classes and if the Sequence/Analysis files are local, uploads
    them to an S3 bucket.
    """

    # Upload local files if input directory is not already an S3 bucket.
    if should_upload_to_s3:
        export_local_files_to_s3(dataset_metadata, input_directory)

    # Upload to DB
    print(f"Committing dataset to Ledger DB")
    session_maker = DBSessionMaker()
    dataset_metadata.export_to_database(session_maker)
    print(f"Completed updating Ledger DB with dataset metadata")


def export_local_files_to_s3(dataset_metadata: OldDatasetMetadata, input_directory):
    """
    Export local files that are directly referenced in the metadata schema to S3.
    """
    print(f"Uploading local files related to the project to S3")

    # Upload FASTQs into S3 and backpopulate dataset_metadata sequence_files
    for filename in listdir(input_directory):
        source_location = os.path.join(input_directory, filename)

        # Just upload the file if this is not a FASTQ file
        if ".json" in filename:
            export_file_to_s3(source_location, filename)
            continue

        # Find a SequenceFile object with similar filename and update it with the S3 URI
        # if it is a FASTQ file
        for sequence_file_id, sequence_file in dataset_metadata.sequence_files.items():
            if sequence_file.filename == filename:
                s3_uri = export_file_to_s3(source_location, filename)
                sequence_file.s3_uri = s3_uri
                dataset_metadata.sequence_files[sequence_file_id] = sequence_file
                continue

        # Find an AnalysisFile object with similar filename and update it with the S3
        # URI if it is not a FASTQ or a JSON file (likely a Zarr or BAM file).
        for analysis_file_id, analysis_file in dataset_metadata.analysis_files.items():
            if analysis_file.filename == filename:
                s3_uri = export_file_to_s3(source_location, filename)
                sequence_file.s3_uri = s3_uri
                dataset_metadata.analysis_files[analysis_file_id] = analysis_file
                continue

    print(f"Completed uploading all files to S3")


def export_file_to_s3(full_source_file_path, filename):
    s3_uri = f"s3://{BUCKET_NAME}/{PREFIX}/{filename}"

    try:
        S3_CLIENT.head_object(Bucket=BUCKET_NAME, Key=PREFIX + filename)
        print(f"Found file {filename} so not re-uploading to S3.")
        return s3_uri
    except Exception:
        pass

    print(f"Uploading file {filename}")
    try:
        S3_CLIENT.upload_file(full_source_file_path, BUCKET_NAME, PREFIX + filename)
    except ClientError as e:
        print(e)
        return False

    return s3_uri


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_directory",
        nargs="+",
        required=False,
        help="A data directory contains all files that were a part of a single DCP 1.0 "
        "project to be transformed into a valid DCP 2.0 metadata and inputted into"
        " the DCP 2.0 ledger.",
    )
    parser.add_argument(
        "-b",
        "--bucket",
        nargs="+",
        required=False,
        help="The bucket into which files should be transfered if transferring local "
        "files as part of the project's metadata migration.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        nargs="+",
        required=False,
        help="The bucket's prefix into which files should be transfered if transferring"
        " local files as part of the project's metadata migration.",
    )

    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
        if "s3" in input_directory:
            # Ensure that input directory ends in a slash
            if not input_directory[-1] == "/":
                print(f"ERROR: Input directory, if an S3 bucket, must end in a slash!")
                sys.exit()

            # Ensure a new bucket and prefix are not given because the metadata will
            # already be entirely read from the S3 bucket. No files will be uploaded to
            # S3.
            if arguments.bucket or arguments.prefix:
                print(
                    f"ERROR: Cannot entire a bucket or prefix will designating that the"
                    f" metadata schema will be read from an S3 bucket."
                )

    if arguments.bucket:
        BUCKET_NAME = arguments.bucket[0]
    if arguments.prefix:
        PREFIX = arguments.prefix[0]

    old_metadata = generate_metadata_structure(input_directory)
    export_old_metadata_to_s3_orm(
        old_metadata, input_directory, "s3" not in input_directory
    )
