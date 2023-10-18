# Extracts CXG file metadata and outputs as json, flattening any dictionary-nested properties into dot-deliminted json key names.
# Can be used for a single CXG file or multiple CXG files using `--url` and `--urls-file`, respectively.
# For usage: `python cxg_inspect.py --help`
import contextlib

# To obtain URLs for CXG files "in the wild":
# - PSQL: select distinct s3_uri from dataset_artifact where filetype= 'CXG' order by created_at;
# - aws s3 ls --recursive s3://cellxgene-data-wrangling-prod/ | grep .cxg | awk '{ print $4 }' | awk -F / '{ print "s3://cellxgene-data-wrangling-prod/"$1"/"$2 }' | sort -u
import json
import os
import pprint
import sys
from typing import Optional

import click

from server.common.config.app_config import AppConfig
from server.common.utils.data_locator import DataLocator
from server.dataset.cxg_dataset import CxgDataset


@click.command()
@click.option("--url", "-u", type=str, help="CXG file URL. E.g. `file:/my_dataset.cxg`, or `s3://fba8204.cxg`.")
@click.option(
    "--urls-file",
    "-f",
    type=str,
    help="File containing a list of CXG file URLs, separated by newlines. " "Output will conform to jsonl format.",
)
@click.option("--limit", "-l", type=int)
def inspect_cxg(url: Optional[str], urls_file: Optional[str], limit=None) -> None:
    urls = []
    if url:
        urls.append(url)
    if urls_file:
        with open(os.path.expanduser(urls_file), "rt") as f:
            urls.extend([f.strip() for f in f.readlines()])

    for cxg_props in filter(None, (_inspect_cxg(url) for url in urls[:limit])):
        print(json.dumps(cxg_props))


def _flatten_dict(d: dict, root_prop: Optional[str] = None):
    if d is None:
        return {}
    flattened_props = {}
    for k, v in d.items():
        fq_key = f"{root_prop}.{k}" if root_prop else k  # fully-qualified key

        # try to decode nested json dict
        if type(v) is str:
            with contextlib.suppress(Exception):
                v = json.loads(v)

        if type(v) is dict:
            flattened_props.update(**_flatten_dict(v, fq_key))
        else:
            flattened_props[fq_key] = v
    return flattened_props


def _inspect_cxg(dataset_url: str) -> Optional[dict]:
    dl = DataLocator(dataset_url)
    if not dl.exists():
        sys.stderr.write(f"NOT FOUND: {dataset_url}\n")
        return None

    try:
        sys.stderr.write(f"Inspecting {dataset_url}\n")
        cxg = CxgDataset(dl, app_config=AppConfig())
        metadata = cxg.open_array("cxg_group_metadata").meta
        return dict(url=dataset_url, **_flatten_dict(metadata))
    except Exception:
        sys.stderr.write(f"No metadata for {dataset_url}\n")
        return dict(url=dataset_url)


# For S3 URLs, these env vars must be set appropriately: AWS_PROFILE, DEPLOYMENT_STAGE, AWS_REGION
if __name__ == "__main__":
    pp = pprint.PrettyPrinter(width=120, sort_dicts=False)
    CxgDataset.set_tiledb_context({"vfs.s3.region": os.getenv("AWS_REGION", "us-west-2")})
    inspect_cxg()
