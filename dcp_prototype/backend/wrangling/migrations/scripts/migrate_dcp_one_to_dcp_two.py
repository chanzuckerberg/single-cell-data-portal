import sys

sys.path.insert(0, "")  # noqa
import boto3
from argparse import ArgumentParser
from urllib.parse import urlparse
import os.path
import json

from dcp_prototype.backend.wrangling.migrations.utils.migration_utils import DatasetMetadata, combine_projects

from dcp_prototype.backend.wrangling.migrations.utils.gather_dcp1_data import (
    generate_metadata_structure_from_targz,
    generate_metadata_structure_from_s3_uri,
    generate_metadata_structure_from_dir,
)


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input-source",
        required=True,
        help="input_source contains all metadata files that were part of a single DCP 1.0 project. It may be one of "
             "the following:\n1. An s3 bucket containing json files.\n2. A directory containing json files.\n3. An s3 "
             "object in the tar.gz format containing the json files.\n4.Or a local tar.gz file containing the json "
             "files.",
    )
    parser.add_argument(
        "-t", "--threads", required=False, default=1, type=int, help="Number of threads to use when reading files"
    )
    parser.add_argument(
        "-o",
        "--output-file",
        default=None,
        help="If provided, the JSON output will be save at this file. If not provided, the JSON output will go to "
             "stdout.",
    )
    parser.add_argument(
        "-a",
        "--append",
        action="store_true",
        help="If provided, the JSON output will be save at this file. If not provided, the JSON output will go to "
             "stdout.",
    )
    arguments = parser.parse_args()

    input_source = arguments.input_source
    output_file = arguments.output_file
    num_threads = arguments.threads

    if "s3" in input_source:
        bucket = urlparse(input_source).netloc
        path = urlparse(input_source).path
        if path.startswith("/"):
            path = path[1:]

        if input_source.endswith(".tar.gz"):
            input_tarfile = os.path.basename(path)
            if os.path.exists(input_tarfile):
                print("ERROR: Cannot overwrite {input_tarfile}")
                sys.exit(1)
            s3_client = boto3.client("s3")
            s3_client.download_file(bucket, path, input_tarfile)
            dataset_metadata = DatasetMetadata()
            generate_metadata_structure_from_targz(input_tarfile, dataset_metadata)
        else:
            if not input_source.endswith("/"):
                input_source += "/"

            dataset_metadata = DatasetMetadata()
            generate_metadata_structure_from_s3_uri(input_source, num_threads, dataset_metadata)

    else:
        dataset_metadata = DatasetMetadata()
        if input_source.endswith(".tar.gz"):
            input_tarfile = input_source
            generate_metadata_structure_from_targz(input_tarfile, dataset_metadata)
        else:
            generate_metadata_structure_from_dir(input_source, num_threads, dataset_metadata)

    dataset_metadata.process()
    result_artifact = dataset_metadata.to_dict()

    okay = dataset_metadata.validate(result_artifact)
    if not okay:
        print("Result failed schema validation")
        sys.exit(1)

    for key, value in dataset_metadata.missing.items():
        print("MISSING:", key, value)

    out_project = result_artifact.get("projects")[0]
    title = out_project.get("title")
    if output_file:
        if arguments.append and os.path.exists(output_file):
            result_artifact = combine_projects(output_file, result_artifact)
        else:
            print(f"Write project {title} to {output_file}")

        okay = dataset_metadata.validate(result_artifact)
        if not okay:
            print("Result failed schema validation")
            sys.exit(1)

        with open(output_file, "w") as json_file:
            json_file.write(json.dumps(result_artifact, indent=2))

    else:
        print(json.dumps(result_artifact, indent=2))


if __name__ == "__main__":
    main()
