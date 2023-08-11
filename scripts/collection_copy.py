#!/usr/bin/env python
import json
import logging
import os
import sys

import click
import requests  # type: ignore

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

logging.basicConfig()
logger = logging.getLogger(__name__)


@click.group()
@click.option("--deployment", default="test", show_default=True, help="The name of the deployment to target.")
@click.pass_context
def cli(ctx, deployment):
    os.environ["DEPLOYMENT_STAGE"] = deployment
    api_base = {
        "dev": "https://api.dev.single-cell.czi.technology/db/v1",
        "staging": "https://api.stage.single-cell.czi.technology/dp/v1",
    }
    ctx.obj["api_base"] = api_base[deployment]


@cli.command()
@click.argument("id")
@click.argument("cookie")
@click.pass_context
def copy_collection(ctx, uuid, cookie):
    """
    Copy a collection from prod to the passed in deployment. Note that the collection owner
    will be set based on the passed in cookie. To retrieve a cookie, sign in (through the website) to the deployment
    you are copying the collection to and copy the cookie (starting with cxguser=) from the request header
    """

    response = requests.get(f"https://api.cellxgene.cziscience.com/dp/v1/collections/{uuid}")
    body = json.loads(response._content)
    dataset_assets = []
    for dataset in body["datasets"]:
        for asset in dataset["dataset_assets"]:
            if asset["filetype"] == "H5AD":
                dataset_assets.append({"asset_id": asset["id"], "dataset_id": asset["dataset_id"]})

    collection_metadata = {
        "contact_email": body["contact_email"],
        "contact_name": body["contact_name"],
        "data_submission_policy_version": body["data_submission_policy_version"],
        "description": body["description"],
        "links": [],
        "name": body["name"],
    }

    for link in body["links"]:
        collection_metadata["links"].append(
            {"link_name": link["link_name"], "link_type": link["link_type"], "link_url": link["link_url"]}
        )

    headers = {"Content-Type": "application/json", "Cookie": cookie, "accept": "application/json"}
    response = requests.post(
        f"{ctx.obj['api_base']}/collections/", headers=headers, data=json.dumps(collection_metadata)
    )

    new_collection_id = json.loads(response._content)["collection_id"]

    for asset in dataset_assets:
        response = requests.post(
            f"https://api.cellxgene.cziscience.com/dp/v1/datasets/{asset['dataset_id']}/asset/{asset['asset_id']}"
        )
        presigned_s3_uri = json.loads(response._content)["presigned_url"]
        dataset_body = {"dataset_id": "", "url": presigned_s3_uri}
        response = requests.post(
            f"{ctx.obj['api_base']}/collections/{new_collection_id}/upload-links",
            headers=headers,
            data=json.dumps(dataset_body),
        )
        dataset_id = json.loads(response._content)
        click.echo(f"New Collection_id: {new_collection_id}, new dataset_id: {dataset_id}")


if __name__ == "__main__":
    cli(obj={})
