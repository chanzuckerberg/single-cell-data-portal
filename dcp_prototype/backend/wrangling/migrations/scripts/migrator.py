import sys

sys.path.insert(0, "")
from argparse import ArgumentParser
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.dataset_metadata import (
    DatasetMetadata,
)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_file",
        nargs=1,
        required=True,
        help="A XLSX file that was part of DCP 1.0 to be transformed into a valid DCP 2.0 spreadsheet.",
    )

    arguments = parser.parse_args()
    input_file = arguments.input_file[0]

    from pandas import read_csv

    metadata_tsv = read_csv(
        input_file, sep="\t", header=0, error_bad_lines=False, encoding="latin-1"
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

    print(
        f"Total rows: {num_rows}. Total filenames: {len(file_names)}. Total unique filenames: {len(set(file_names))}"
    )
    print(
        f"Total rows: {num_rows}. Total file_uuids: {len(file_uuids)}. Total unique file_uuids: {len(set(file_uuids))}"
    )
    dataset_metadata.export_to_spreadsheet()
