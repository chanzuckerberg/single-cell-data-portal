import sys

sys.path.insert(0, "")  # noqa
from argparse import ArgumentParser
from pandas.io.json import json_normalize
from os import listdir
import os.path
import json
from flatten_json import flatten


def generate_metadata_tsv_name(project_name):
    """
    Given a project name, generates a name for the spreadsheet that will be generated
    based on this project's metadata.
    """
    replace_spaces = project_name.replace(" ", "_")
    replace_slashes = replace_spaces.replace("/", "_")
    return replace_slashes + "_metadata.xlsx"


def generate_single_tsv(input_filename):
    """
    For a given filename which should be a TSV file containing all of the flattened
    metadata representing a single project from DCP 1.0, transforms that project's
    metadata into the DCP 2.0 metadata schema and projects it into an XLSX spreadsheet
    for manual validation.
    """

    dataset_metadata = DatasetMetadata()

    with open(input_filename) as json_file_object:
        input_file_contents = json_file_object.read()

    tsv_generated_data_frame = json_normalize(
        flatten(json.loads(input_file_contents), separator=".")
    )
    tsv_generated_data_frame.to_excel(
        generate_metadata_tsv_name(input_filename), index=False
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        nargs=1,
        required=False,
        help="A TSV file that was part of DCP 1.0 to be transformed into a valid DCP "
        "2.0 spreadsheet.",
    )
    parser.add_argument(
        "-d",
        "--input_directory",
        nargs=1,
        required=False,
        help="A data directory contains TSV files that were a part of DCP 1.0 to be "
        "transformed into a valid DCP 2.0 spreadsheets.",
    )

    arguments = parser.parse_args()

    if arguments.input_file:
        input_file = arguments.input_file[0]
        generate_single_tsv(input_file)

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
        print(f"Files in directory to parse: {listdir(input_directory)}")
        for file in listdir(input_directory):
            full_file_path = os.path.join(input_directory, file)
            if os.path.isfile(full_file_path):
                print(f"Working with input spreadsheet at {full_file_path}")
                generate_single_ts(full_file_path)
