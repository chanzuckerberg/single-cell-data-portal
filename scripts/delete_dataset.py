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
s3 = boto3.resource("s3")


@click.command()
@click.argument("uuid")
@click.option(
    "--deployment", default=os.getenv("DEPLOYMENT_STAGE", "test"), help="The name of the deployment to target."
)
def delete_dataset(uuid, deployment):
    """Delete a dataset from Cellxgene. You must first SSH into the target deployment using `make db/tunnel` before
    running."""
    DBSessionMaker(get_database_uri())
    with db_session_manager() as session:
        dataset = Dataset.get(session, uuid, include_tombstones=True)
        if dataset is not None:
            print(json.dumps(dataset.to_dict(), sort_keys=True, indent=2, cls=CustomJSONEncoder))
            delete_deployment_directories(dataset.deployment_directories, deployment)
            dataset.dataset_and_asset_deletion()
            dataset.delete()
            dataset = Dataset.get(session, uuid, include_tombstones=True)
        if dataset is None:
            click.echo(f"Deleted: {uuid}")
            exit(0)
        else:
            click.echo(f"Failed to delete: {uuid}")
            exit(1)


def get_database_uri() -> str:
    uri = urlparse(CorporaDbConfig().database_uri)
    uri = uri._replace(netloc="@".join([uri[1].split("@")[0], "localhost:5432"]))
    return uri.geturl()


def delete_deployment_directories(deployment_directories, deployment):
    s3 = boto3.resource("s3")
    bucket_name = f"hosted-cellxgene-{deployment}"
    bucket = s3.Bucket(bucket_name)
    for deployment_directory in deployment_directories:
        object_name = urlparse(deployment_directory.url).path.split("/", 2)[2]
        logger.info(f"Deleting all files in bucket {bucket_name} under {object_name}.")
        bucket.objects.filter(Prefix=object_name).delete()


if __name__ == "__main__":
    delete_dataset()
