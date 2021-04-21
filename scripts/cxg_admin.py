#!/usr/bin/env python
import boto3
import click
import json
import logging
import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.corpora_config import CorporaDbConfig
from backend.corpora.common.utils.json import CustomJSONEncoder
from backend.corpora.common.utils.db_session import db_session_manager, DBSessionMaker
from backend.corpora.common.entities.dataset import Dataset

from urllib.parse import urlparse

logging.basicConfig()
logger = logging.getLogger(__name__)

os.environ["CORPORA_LOCAL_DEV"] = "1"


@click.group()
@click.option("--deployment", default="test", show_default=True, help="The name of the deployment to target.")
@click.pass_context
def cli(ctx, deployment):
    os.environ["DEPLOYMENT_STAGE"] = deployment
    ctx.obj["deployment"] = deployment
    DBSessionMaker(get_database_uri())


@cli.command()
@click.argument("uuid")
@click.pass_context
def delete_dataset(ctx, uuid):
    """Delete a dataset from Cellxgene. You must first SSH into the target deployment using `make db/tunnel` before
    running."""

    with db_session_manager() as session:
        dataset = Dataset.get(session, uuid, include_tombstones=True)
        if dataset is not None:
            print(
                json.dumps(dataset.to_dict(remove_attr=["collection"]), sort_keys=True, indent=2, cls=CustomJSONEncoder)
            )
            click.confirm(
                f"Are you sure you want to delete the dataset:{uuid} from cellxgene:{ctx.obj['deployment']}?",
                abort=True,
            )
            delete_deployment_directories(dataset.deployment_directories, ctx.obj["deployment"])
            dataset.dataset_and_asset_deletion()
            dataset.delete()
            dataset = Dataset.get(session, uuid, include_tombstones=True)
            if dataset is None:
                click.echo(f"Deleted dataset:{uuid}")
                exit(0)
            else:
                click.echo(f"Failed to delete dataset:{uuid}")
                exit(1)
        else:
            click.echo(f"Dataset:{uuid} not found!")
            exit(0)


def get_database_uri() -> str:
    uri = urlparse(CorporaDbConfig().database_uri)
    uri = uri._replace(netloc="@".join([uri[1].split("@")[0], "localhost:5432"]))
    return uri.geturl()


def delete_deployment_directories(deployment_directories, deployment):
    s3 = boto3.resource("s3", endpoint_url=os.getenv("BOTO_ENDPOINT_URL"))
    bucket_name = f"hosted-cellxgene-{deployment}"
    bucket = s3.Bucket(bucket_name)
    for deployment_directory in deployment_directories:
        object_name = urlparse(deployment_directory.url).path.split("/", 2)[2]
        logger.info(f"Deleting all files in bucket {bucket_name} under {object_name}.")
        bucket.objects.filter(Prefix=object_name).delete()


if __name__ == "__main__":
    cli(obj={})
