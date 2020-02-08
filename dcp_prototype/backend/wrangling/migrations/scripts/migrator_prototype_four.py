import sys

sys.path.insert(0, "")  # noqa
from argparse import ArgumentParser
from pandas.io.json import json_normalize
from os import listdir, mkdir
import os.path
import json
from flatten_json import flatten
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.dataset_metadata import (  # noqa
    DatasetMetadata,
)
import hashlib
from shutil import copyfile, rmtree

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


def generate_tsv_from_bundle(input_directory):
    """
    For a given filename which should be a TSV file containing all of the flattened
    metadata representing a single project from DCP 1.0, transforms that project's
    metadata into the DCP 2.0 metadata schema and projects it into an XLSX spreadsheet
    for manual validation.
    """

    dataset_metadata = DatasetMetadata()

    ordered_files = order_file_list(listdir(input_directory))

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


def squish_files(directory, new_directory_path, checksums=[], clear_if_exists=True):

    if os.path.isdir(new_directory_path) and clear_if_exists:
        rmtree(new_directory_path)
    elif not os.path.isdir(new_directory_path):
        mkdir(new_directory_path)

    current_path = directory
    for file in listdir(directory):
        source_location = os.path.join(current_path, file)
        if file.startswith("."):
            continue
        if os.path.isfile(source_location):
            dest_location = os.path.join(new_directory_path, file)
            filename = file

            with open(source_location, "rb") as file_reader:
                contents = file_reader.read()

            checksum = hashlib.md5(contents).hexdigest()
            print(f"Processing file: {source_location} with checksum {checksum}")
            if checksum not in checksums:
                checksums.append(checksum)

                if ".json" in file:
                    for file_type in COUNT_BY_FILE_TYPE.keys():
                        if file_type in file:
                            if file_type == "donor_organism":
                                print(f"NUMBER: {COUNT_BY_FILE_TYPE[file_type]}")
                            filename = (
                                f"{file_type}_{str(COUNT_BY_FILE_TYPE[file_type])}.json"
                            )
                            COUNT_BY_FILE_TYPE[file_type] += 1
                            if file_type == "donor_organism":
                                print(filename)


                copyfile(source_location, os.path.join(new_directory_path, filename))
            else:
                print(
                    f"Skipping over file {file} because it has already been copied over."
                )
        elif os.path.isdir(source_location):
            print(f"Processing directory {os.path.join(current_path, file)}")
            squish_files(
                os.path.join(current_path, file), new_directory_path, checksums, False
            )


def should_process_file(full_file_path):
    if "analysis" in full_file_path:
        return False
    if "zarr" in full_file_path:
        return False
    if ".json" not in full_file_path:
        return False
    return True


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-d",
        "--input_directory",
        nargs="+",
        required=False,
        help="A data directory contains TSV files that were a part of DCP 1.0 to be "
        "transformed into a valid DCP 2.0 spreadsheets.",
    )
    parser.add_argument(
        "-o", "--output_directory", nargs="+", required=False,
    )

    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
    if arguments.output_directory:
        output_directory = arguments.output_directory[0]
    print(f"Files in directory to parse: {listdir(input_directory)}")
    #squish_files(input_directory, output_directory)
    generate_tsv_from_bundle(input_directory)
