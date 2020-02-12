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
    dataset_metadata = OldDatasetMetadata()

    BUCKET_NAME = urlparse(s3_uri).netloc
    PREFIX = urlparse(s3_uri).path
    if PREFIX.startswith("/"):
        PREFIX = PREFIX[1:]
    response = S3_CLIENT.list_objects(Bucket=BUCKET_NAME, Prefix=PREFIX)

    bucket_objects = response.get("Contents")
    file_list = []
    for object in bucket_objects:
        object_filename = object.get("Key")
        file_list.append(object_filename.split("/")[-1])
    file_list = order_file_list(file_list)

    print(f"Files in directory to parse: {file_list}")

    for filename in file_list:
        if not should_process_file(filename):
            continue

        file_prefix = PREFIX + "/" + filename

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

    dataset_metadata.save()

    print(f"Processing dataset: {dataset_metadata.project_name}")
    dataset_metadata.export_to_spreadsheet(
        generate_metadata_tsv_name(dataset_metadata.project_name)
    )
    print(
        f"Completed generating spreadsheet: "
        f"{generate_metadata_tsv_name(dataset_metadata.project_name)}"
    )

    return dataset_metadata


def generate_metadata_structure(input_directory):
    """
    For a given filename which should be a TSV file containing all of the flattened
    metadata representing a single project from DCP 1.0, transforms that project's
    metadata into the DCP 2.0 metadata schema and projects it into an XLSX spreadsheet
    for manual validation.
    """

    if "s3" in input_directory:
        return generate_metadata_structure_from_s3_uri(input_directory)

    dataset_metadata = OldDatasetMetadata()

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

    print(f"Processing dataset: {dataset_metadata.project_name}")
    dataset_metadata.export_to_spreadsheet(
        generate_metadata_tsv_name(dataset_metadata.project_name)
    )
    print(
        f"Completed generating spreadsheet: "
        f"{generate_metadata_tsv_name(dataset_metadata.project_name)}"
    )

    return dataset_metadata


def should_process_file(full_file_path):
    if "analysis" in full_file_path:
        return False
    if "zarr" in full_file_path:
        return False
    if ".json" not in full_file_path:
        return False
    return True


def export_old_metadata_to_s3_orm(
    dataset_metadata: OldDatasetMetadata, input_directory
):

    # Upload FASTQs into S3 and backpopulate dataset_metadata sequence_files
    for filename in listdir(input_directory):
        source_location = os.path.join(input_directory, filename)

        # Just upload the file if this is not a FASTQ file
        if ".fastq" not in filename:
            export_file_to_s3(source_location, filename)
            continue

        # Find SequenceFile object with similar filename and update it with the S3 URI
        # if it is a FASTQ file
        for sequence_file_id, sequence_file in dataset_metadata.old_files.items():
            if sequence_file.filename == filename:
                s3_uri = export_file_to_s3(source_location, filename)
                sequence_file.set_s3_uri(s3_uri)
                dataset_metadata.old_files[sequence_file_id] = sequence_file
                continue
    print(f"Completed uploading all files to S3")

    dataset_metadata.convert_to_new_entities()

    # Upload to DB
    print(f"Commit dataset to DB")
    session_maker = DBSessionMaker()
    session = session_maker.session()

    for biosample_prep in dataset_metadata.biosample_preps:
        session.add(biosample_prep)
    for sequence_protocol in dataset_metadata.sequencing_protocols:
        session.add(sequence_protocol)
    for project in dataset_metadata.projects:
        session.add(project)
    for contributor in dataset_metadata.contributors:
        session.add(contributor)
    for sequence_file in dataset_metadata.sequence_files:
        session.add(sequence_file)
    for library_prep in dataset_metadata.library_preps:
        session.add(library_prep)
    for library in dataset_metadata.libraries:
        session.add(library)

    session.commit()


def export_file_to_s3(full_source_file_path, filename):
    s3_uri = f"s3://{BUCKET_NAME}/{filename}"

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

    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]

    old_metadata = generate_metadata_structure(input_directory)
    export_old_metadata_to_s3_orm(old_metadata, input_directory)
