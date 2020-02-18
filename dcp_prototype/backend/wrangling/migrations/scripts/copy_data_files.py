# import concurrent
import json
import os

import boto3
from argparse import ArgumentParser


def copy_between_s3_buckets(file_name, project):
    s3 = boto3.resource('s3')
    copy_source = {
        'Bucket': 'org-hca-dss-prod',
        'Key': f'files/{file_name}'
    }
    s3.meta.client.copy(copy_source, 'dunitz-prod-copy', f'{project}/{file_name}')


def list_data_files_for_project(input_directory, project_name):
    project_data_files = []
    target_directory_path = os.path.join(input_directory, project_name, 'bundle_manifests')
    for file in os.listdir(target_directory_path):
        # skip hidden files
        if file.startswith("."):
            continue
        file_path = os.path.join(target_directory_path, file)
        print(file_path)
        with open(file_path, "rb") as file_reader:
            contents = file_reader.read()
            files = json.loads(contents)
            for file_info in files['files']:
                if file_info["content-type"].find("'application/json; dcp-type=\"metadata*\"'")
    return []


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_directory",
        nargs="+",
        required=True,
        help="A data directory containing bundles with files that contain both the "
        "metadata and the data of a single project.",
    )
    arguments = parser.parse_args()

    if arguments.input_directory:
        input_directory = arguments.input_directory[0]
    for project in os.listdir(input_directory):
        project_file_list = list_data_files_for_project(input_directory, project)
    # dispatch_executor_class = concurrent.futures.ThreadPoolExecutor
    projects = []
    # with dispatch_executor_class(max_workers=50) as executor:
    #     for project in projects:
    #         futures = []
    #         f = executor.submit(list_data_files_for_project, project)
    #         futures.append(f)
    #
    #         for future in concurrent.futures.as_completed(futures):
    #             try:
    #                 extract_result = future.result()
    #                 for file in extract_result:
    #                     file_name = f"{file['uuid']}.{file['version']}"
    #                     copy_between_s3_buckets(file_name, project)
    #             except Exception as e:
    #                 print(f"Something went wrong: {e}")
