import sys

sys.path.insert(0, "")
from argparse import ArgumentParser
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.dataset_metadata import (
    DatasetMetadata,
)
from pandas import read_csv
from os import listdir
import os.path


def generate_metadata_spreadsheet_name(project_name):
    replace_spaces = project_name.replace(" ", "_")
    replace_slashes = replace_spaces.replace("/", "_")
    return replace_slashes + "_metadata.xlsx"


def generate_single_spreadsheet(input_filename):
    metadata_tsv = read_csv(
        input_filename, sep="\t", header=0, error_bad_lines=False, encoding="latin-1"
    )
    num_rows = 0
    file_names = []
    file_uuids = []
    dataset_metadata = DatasetMetadata()
    for row in metadata_tsv.iterrows():
        metadata_row = row[1]
        num_rows += 1
        if metadata_row.get("sequence_file.file_core.format") == "fastq.gz":
            file_names.append(metadata_row.get("*.file_core.file_name"))
            file_uuids.append(metadata_row.get("file_uuid"))
            dataset_metadata.parse_row_of_metadata(metadata_row)
    dataset_metadata.save()
    print(f"Processing dataset: {dataset_metadata.project_name}")
    print(
        f"Total rows: {num_rows}. Total filenames: {len(file_names)}. Total unique filenames: {len(set(file_names))}"
    )
    print(
        f"Total rows: {num_rows}. Total file_uuids: {len(file_uuids)}. Total unique file_uuids: {len(set(file_uuids))}"
    )
    dataset_metadata.export_to_spreadsheet(
        generate_metadata_spreadsheet_name(dataset_metadata.project_name)
    )
    print(
        f"Completed generating spreadsheet: {generate_metadata_spreadsheet_name(dataset_metadata.project_name)}"
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        nargs=1,
        required=False,
        help="A XLSX file that was part of DCP 1.0 to be transformed into a valid DCP 2.0 spreadsheet.",
    )
    parser.add_argument(
        "-d",
        "--input_directory",
        nargs=1,
        required=False,
        help="A data directory contains XLSX files that were a part of DCP 1.0 to be transformed into a valid DCP 2.0 spreadsheets.",
    )

    arguments = parser.parse_args()

    if arguments.input_file:
        input_file = arguments.input_file[0]
        generate_single_spreadsheet(input_file)

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
        print(f"Files in directory to parse: {listdir(input_directory)}")
        for file in listdir(input_directory):
            full_file_path = os.path.join(input_directory, file)
            if os.path.isfile(full_file_path):
                print(f"Working with input spreadsheet at {full_file_path}")
                generate_single_spreadsheet(full_file_path)

