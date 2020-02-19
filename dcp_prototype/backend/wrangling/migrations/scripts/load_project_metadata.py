import json
from argparse import ArgumentParser
import sys

sys.path.insert(0, "")  # noqa
from dcp_prototype.backend.wrangling.migrations.scripts.squish_files import squish_files, reset_counts


def load_project_metadata(input_directory, output_directory, count_file, project_file):
    with open(project_file) as json_file:
        data = json.load(json_file)
        for project in data:
            project = project.replace(" ", "-").replace("/", "-").replace(".", "")
            full_input_directory = f"{input_directory}/{project}/bundles"
            full_output_directory = f"{output_directory}/{project}"
            reset_counts(count_file)
            squish_files(full_input_directory, full_output_directory, count_file)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_directory",
        nargs="+",
        required=True,
        help="A data directory containing the metadata files for one or more projects"
    )
    parser.add_argument(
        "-o",
        "--output_directory",
        nargs="+",
        required=True,
        help="An output directory to which the files from the input directory will be "
             "copied to and de-duped. File names will be retained except for the suffix"
             " that is a numbering of the file.",
    )
    parser.add_argument(
        "-c",
        "--count_file",
        nargs="+",
        required=True,
        help="A count file that keeps track of the number of a particular file that the"
             " script has processed. This file allows for multiple runs of the script "
             "while retaining previous counts.",
    )
    parser.add_argument(
        "-p",
        "--project_file",
        nargs="+",
        required=True,
        help="A file containing a list of projects to load metadata for",
    )

    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
    if arguments.output_directory:
        output_directory = arguments.output_directory[0]
    if arguments.count_file:
        count_file = arguments.count_file[0]
    if arguments.project_file:
        project_file = arguments.project_file[0]

    load_project_metadata(input_directory, output_directory, count_file, project_file)
