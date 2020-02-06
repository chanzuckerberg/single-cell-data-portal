import sys

sys.path.insert(0, "")  # noqa
from argparse import ArgumentParser
from pandas.io.json import json_normalize
from os import listdir
import os.path
import json
from flatten_json import flatten
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.dataset_metadata import (  # noqa
    DatasetMetadata,
)

MAP_FROM_FILE_NAME_TO_ENTITY_TYPE = {
    "donor_organism": "biosample_prep",
    "library_preparation_protocol": "library_preparation_protocol",
    "project": "project",
    "sequence_file": "sequence_file",
    "sequencing_protocol": "sequencing_protocol",
    "specimen_from_organism" : "specimen",
    "cell_suspension" : "other_biosample",
    "links": "links",
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
    for key in MAP_FROM_FILE_NAME_TO_ENTITY_TYPE.keys():
        if key in filename:
            return MAP_FROM_FILE_NAME_TO_ENTITY_TYPE[key]

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
        if file not in ordered_file_list:
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

    links_json_full_file_path = os.path.join(input_directory, "links.json")
    print(f"Working with JSON input file at {links_json_full_file_path}")
    with open(links_json_full_file_path) as json_file_object:
        input_file_contents = json_file_object.read()
    tsv_generated_data_frame = json_normalize(
        flatten(json.loads(input_file_contents), separator=".")
    )
    for row in tsv_generated_data_frame.iterrows():
        metadata_row = row[1]
        dataset_metadata.parse_flattened_row_of_json(
            metadata_row, get_entity_type(links_json_full_file_path)
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


def should_process_file(full_file_path):
    if "analysis" in full_file_path:
        return False
    if "zarr" in full_file_path:
        return False
    if ".json" not in full_file_path:
        return False
    if "links.json" in full_file_path:
        return False
    return True


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-d",
        "--input_directory",
        nargs=1,
        required=False,
        help="A data directory contains TSV files that were a part of DCP 1.0 to be "
        "transformed into a valid DCP 2.0 spreadsheets.",
    )

    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
        print(f"Files in directory to parse: {listdir(input_directory)}")
        generate_tsv_from_bundle(input_directory)
