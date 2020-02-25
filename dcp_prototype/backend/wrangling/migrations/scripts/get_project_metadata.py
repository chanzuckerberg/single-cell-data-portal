import json
from argparse import ArgumentParser
from dcplib.etl import DSSExtractor


def get_project_metadata(project_file, max_workers):
    """
    Get metadata (.json files) for all projects listed in project_file. Storing the data in ./hold_data/{project_name}/
    :param project_file: json file contanining a list of projects
    :param max_workers: Max numbers of workers for the extractor
    """
    with open(project_file) as json_file:
        data = json.load(json_file)
        for project in data:
            query = {
                "query": {
                    "bool": {"must": [{"term": {"files.project_json.project_core.project_short_name": f"{project}"}}]}
                }
            }
            project = project.replace(" ", "-").replace("/", "-").replace(".", "")
            print(project)
            DSSExtractor(staging_directory=f"./hold_data/{project}").extract(query=query, max_workers=max_workers)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-t", "--threads", nargs="+", required=False, help="The number of threads for extraction, default is 1",
    )
    parser.add_argument(
        "-p",
        "--project_file",
        nargs="+",
        required=True,
        help="A file containing a list of projects to retrieve metadata for",
    )
    arguments = parser.parse_args()
    if arguments.threads:
        threads = int(arguments.threads[0])
    else:
        threads = 1
    if arguments.project_file:
        project_file = arguments.project_file[0]
    get_project_metadata(project_file, threads)
