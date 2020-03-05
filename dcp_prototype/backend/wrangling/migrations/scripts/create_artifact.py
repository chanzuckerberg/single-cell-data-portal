import sys

sys.path.insert(0, "")  # noqa
import argparse
from urllib.parse import urlparse
import boto3
import os
import tempfile
import json

from dcp_prototype.backend.wrangling.migrations.utils.migration_utils import DatasetMetadata, combine_projects

from dcp_prototype.backend.wrangling.migrations.utils.gather_dcp1_data import generate_metadata_structure_from_targz


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_directory",
        help="An s3 prefix contains metadata files that were a part of a single DCP 1.0",
        default="s3://hca-dcp-one-backup-data/dcp-one-copy/",
    )
    parser.add_argument("-o", "--output-file", help="the json output will be save at this file.", required=True)
    arguments = parser.parse_args()

    input_directory = arguments.input_directory
    if input_directory.startswith("s3:"):
        projects = []
        bucket_name = urlparse(input_directory).netloc
        prefix = urlparse(input_directory).path
        if prefix.startswith("/"):
            prefix = prefix[1:]
        client = boto3.client("s3")
        resp = client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")
        for obj in resp["Contents"]:
            key = obj["Key"]
            if key.endswith(".tar.gz"):
                projects.append(os.path.join(key))

        with tempfile.TemporaryDirectory() as dirname:
            for s3name in projects:
                print(f"Processing {s3name}")
                tarfile = os.path.join(dirname, os.path.basename(s3name))
                client.download_file(bucket_name, s3name, tarfile)
                dataset_metadata = DatasetMetadata()
                generate_metadata_structure_from_targz(tarfile, 1, dataset_metadata)
                dataset_metadata.process()
                result_artifact = dataset_metadata.to_dict()

                okay = dataset_metadata.validate(result_artifact)
                if not okay:
                    print("result failed schema validation")
                    sys.exit(1)

                result_artifact = combine_projects(arguments.output_file, result_artifact)
                okay = dataset_metadata.validate(result_artifact)
                if not okay:
                    print("result failed schema validation")
                    sys.exit(1)

                with open(arguments.output_file, "w") as json_file:
                    json_file.write(json.dumps(result_artifact, indent=2))


if __name__ == "__main__":
    main()
