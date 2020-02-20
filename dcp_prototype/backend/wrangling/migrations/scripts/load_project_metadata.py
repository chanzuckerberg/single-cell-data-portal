import concurrent.futures
import json
from argparse import ArgumentParser
import sys

sys.path.insert(0, "")  # noqa
from dcp_prototype.backend.wrangling.migrations.scripts.squish_files import squish_files, reset_counts


def load_project_metadata(input_directory, output_directory, project_file, max_workers, should_reset):
    with open(project_file) as json_file:
        data = json.load(json_file)
        dispatch_executor_class = concurrent.futures.ThreadPoolExecutor
        with dispatch_executor_class(max_workers=max_workers) as executor:
            futures = []
            projects_loaded = 0
            for project in data:
                project = project.replace(" ", "-").replace("/", "-").replace(".", "")
                full_input_directory = f"{input_directory}/{project}/bundles"
                full_output_directory = f"{output_directory}/{project}"
                count_file = f"{input_directory}/{project}/count_file.json"
                if should_reset:
                    reset_counts(count_file)
                f = executor.submit(squish_files, full_input_directory, full_output_directory, count_file,
                                    clear_if_exists=should_reset)
                futures.append(f)

            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                    projects_loaded += 1
                    print(f"{projects_loaded} projects loaded so far (out of {len(data)})")
                except Exception as e:
                    print(f"Something went wrong: {e}")


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
        "-m",
        "--max_workers",
        nargs="+",
        required=False,
        help="The number of threads for extraction, default is 1",
    )
    parser.add_argument(
        "-p",
        "--project_file",
        nargs="+",
        required=True,
        help="A file containing a list of projects to load metadata for",
    )
    parser.add_argument("-r", "--reset", action="store_true")

    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
    if arguments.output_directory:
        output_directory = arguments.output_directory[0]
    if arguments.max_workers:
        max_workers = int(arguments.max_workers[0])
    else:
        max_workers = 1
    if arguments.project_file:
        project_file = arguments.project_file[0]
    should_reset = arguments.reset

    load_project_metadata(input_directory, output_directory, project_file, max_workers, should_reset)
