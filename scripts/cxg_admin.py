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
from backend.corpora.common.corpora_orm import CollectionVisibility, DbCollection
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.entities.collection import Collection

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
            dataset.asset_deletion()
            dataset.deployment_directory_deletion()
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


@cli.command()
@click.argument("collection_uuid")
@click.argument("new_owner")
@click.pass_context
def update_collection_owner(ctx, collection_uuid, new_owner):
    """Update the owner of a cellxgene collection. You must first SSH into the target deployment using
    `make db/tunnel` before running."""

    with db_session_manager() as session:

        key = (collection_uuid, CollectionVisibility.PUBLIC.name)
        collection = Collection.get(session, key)
        collection_name = collection.to_dict()["name"]

        if collection is not None:
            click.confirm(
                f"Are you sure you want to update the owner of the collection:{collection_name}?",
                abort=True,
            )
            collection.update(owner=new_owner)
            collection = Collection.get(session, key)
            if collection and (collection.owner == new_owner):
                click.echo(f"Updated owner of collection:{collection_uuid}. Owner is now {collection.owner}")
                exit(0)
            else:
                click.echo(f"Failed to update owner for collection_uuid:{collection_uuid}")
                exit(0)


@cli.command()
@click.argument("curr_owner")
@click.argument("new_owner")
@click.pass_context
def transfer_collections(ctx, curr_owner, new_owner):
    """Transfer all collections owned bu the curr_owner to the new_owner. You must first SSH into the target
    deployment using `make db/tunnel` before running."""

    with db_session_manager() as session:
        collections = session.query(DbCollection).filter(DbCollection.owner == curr_owner).all()

        if collections is not None:
            click.confirm(
                f"Are you sure you want to update the owner of {len(collections)} collection{'s' if len(collections)>1 else ''} from {curr_owner} to "
                f"{new_owner}?",
                abort=True,
            )
            updated = (
                session.query(DbCollection)
                .filter(DbCollection.owner == curr_owner)
                .update({DbCollection.owner: new_owner})
            )
            if updated > 0:
                click.echo(f"Updated owner of collection for {updated} collections. {new_owner} is now the owner")
                exit(0)
            else:
                click.echo(
                    f"Failed to update owner for collections. {curr_owner} is still the owner of {len(collections)} "
                    f"collections"
                )
                exit(0)


def get_database_uri() -> str:
    uri = urlparse(CorporaDbConfig().database_uri)
    uri = uri._replace(netloc="@".join([uri[1].split("@")[0], "localhost:5432"]))
    return uri.geturl()


if __name__ == "__main__":
    cli(obj={})
