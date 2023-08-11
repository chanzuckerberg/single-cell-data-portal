# Authored by Bruce Martin
# The script is used for checking if our /curation/* endpoints stay performant
# after being hit concurrently in a short period of time
# Context: https://czi-sci.slack.com/archives/C023Q1APASK/p1667438044981199?thread_ts=1667431712.857239&cid=C023Q1APASK
# Note: The Slack link will expire around 11/2/2023 or whatever CZI's Slack retention policy entails

import argparse
import concurrent.futures
import logging
import multiprocessing
import sys
import urllib.parse

import requests  # type: ignore

CXG_BASE_URI = "https://api.cellxgene.cziscience.com/"  # default - can override on the CLI


def load_datasets_manifest_from_CxG(args: argparse.Namespace) -> dict[str, dict]:
    logging.info("Loading datasets manifest from CELLxGENE data portal...")

    # Load all collections and extract dataset_id
    collections = fetch_json(urllib.parse.urljoin(args.api_url, "curation/v1/collections"))
    datasets = {
        dataset["id"]: {
            "collection_id": collection["id"],
            "collection_name": collection["name"],
            "dataset_title": dataset.get("title", ""),  # title is optional in schema
            "dataset_id": dataset["id"],
        }
        for collection in collections
        for dataset in collection["datasets"]
        if dataset["tombstone"] is False  # ignore anything that has been deleted
    }
    logging.info(f"Found {len(datasets)} datasets, in {len(collections)} collections")

    # load per-dataset schema version
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as tp:
        dataset_metadata = tp.map(
            lambda d: fetch_json(
                urllib.parse.urljoin(
                    args.api_url,
                    f"curation/v1/collections/{d['collection_id']}/datasets/{d['dataset_id']}",
                )
            ),
            datasets.values(),
        )
    for d in dataset_metadata:
        datasets[d["id"]].update(
            {
                "schema_version": d["schema_version"],
                "dataset_title": d["title"],
            }
        )

    # Grab the asset URI for each dataset
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as tp:
        dataset_assets = tp.map(
            lambda d: (
                d["dataset_id"],
                fetch_json(
                    urllib.parse.urljoin(
                        args.api_url,
                        f"curation/v1/collections/{d['collection_id']}/datasets/{d['dataset_id']}/assets",
                    )
                ),
            ),
            datasets.values(),
        )

    no_asset_found = []
    for dataset_id, assets in dataset_assets:
        assets_h5ad = [a for a in assets if a["filetype"] == "H5AD"]
        if len(assets_h5ad) == 0:
            logging.error(f"Unable to find H5AD asset for dataset id {dataset_id} - ignoring this dataset")
            no_asset_found.append(dataset_id)
        else:
            asset = assets_h5ad[0]
            datasets[dataset_id].update(
                {
                    "corpora_asset_h5ad_uri": asset["presigned_url"],
                    "asset_h5ad_filesize": asset["filesize"],
                }
            )

    return datasets


def main():
    parser = create_args_parser()
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose > 1 else logging.INFO if args.verbose == 1 else logging.WARNING
    logging.basicConfig(
        format="%(asctime)s %(process)-7s %(levelname)-8s %(message)s",
        level=level,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # code elsewhere assumes base URI ends with a slash
    args.api_url = args.api_url + "/" if not args.api_url.endswith("/") else args.api_url

    while True:
        load_datasets_manifest_from_CxG(args)
        if not args.repeat:
            break

    return 0


def create_args_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="api_thrasher")
    parser.add_argument(
        "--api-url", type=str, required=False, default=CXG_BASE_URI, help=f"Base API URL (default: {CXG_BASE_URI})"
    )
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase logging verbosity")
    parser.add_argument("--repeat", action=argparse.BooleanOptionalAction, default=False, help="Loop until killed")
    parser.add_argument("--max-workers", type=int, default=None, required=False, help="Thread pool size")
    return parser


def fetch_json(url: str) -> dict:
    logging.info(f"fetch {url}")
    response = requests.get(url)
    response.raise_for_status()
    return response.json()


if __name__ == "__main__":
    if multiprocessing.get_start_method(True) != "spawn":
        multiprocessing.set_start_method("spawn", True)

    sys.exit(main())
