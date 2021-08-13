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
from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
    DbCollection,
    DbDataset,
    DatasetArtifactFileType,
    DatasetArtifactType,
    DbDatasetArtifact,
)
from backend.corpora.common.entities import DatasetAsset
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.entities.collection import Collection
from backend.corpora.common.utils.s3_buckets import cxg_bucket

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
@click.argument("uuid")
@click.pass_context
def tombstone_dataset(ctx, uuid):
    """
    Remove a dataset from Cellxgene. This will delete its artifacts/genesets and mark the dataset as tombstoned so
     it no longer shows up in the data portal.
     You must first SSH into the target deployment using `make db/tunnel` before running.
      ./scripts/cxg_admin.py --deployment staging tombstone-dataset "57cf1b53-af10-49e5-9a86-4bc70d0c92b6"

     """

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
            dataset.tombstone_dataset_and_delete_child_objects()
            dataset = Dataset.get(session, uuid)
            tombstoned = Dataset.get(session, uuid, include_tombstones=True)
            if tombstoned and dataset is None:
                click.echo(f"Tombstoned dataset:{uuid}")
                exit(0)
            else:
                click.echo(f"Failed to tombstone dataset:{uuid}")
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
    `make db/tunnel` before running.
    To run (from repo root)
    ./scripts/cxg_admin.py --deployment prod update-collection-owner "$COLLECTION_ID $NEW_OWNER_ID
    """

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
    """Transfer all collections owned by the curr_owner to the new_owner. You must first SSH into the target
    deployment using `make db/tunnel` before running.
    Retrieve user ids from auth0 before running or ping an engineer on the team to check the id of the owner in the database
    To run (from repo root)
    ./scripts/cxg_admin.py --deployment prod transfer-collections $CURR_OWNER_ID $NEW_OWNER_ID
    """

    with db_session_manager() as session:
        collections = session.query(DbCollection).filter(DbCollection.owner == curr_owner).all()
        new_owner_collections_count = len(session.query(DbCollection).filter(DbCollection.owner == new_owner).all())
        if len(collections):
            click.confirm(
                f"Are you sure you want to update the owner of {len(collections)} collection{'s' if len(collections) > 1 else ''} from {curr_owner} to "
                f"{new_owner}?",
                abort=True,
            )
            updated = (
                session.query(DbCollection)
                .filter(DbCollection.owner == curr_owner)
                .update({DbCollection.owner: new_owner})
            )
            session.commit()
            if updated > 0:
                collections = session.query(DbCollection).filter(DbCollection.owner == new_owner).all()
                click.echo(
                    f"{new_owner} previously owned {new_owner_collections_count}, they now own {len(collections)}"
                )
                click.echo(
                    f"Updated owner of collection for {updated} collections. {new_owner} is now the owner of {[[x.name, x.id] for x in collections]}"
                )
                exit(0)
            else:
                click.echo(
                    f"Failed to update owner for collections. {curr_owner} is still the owner of {len(collections)} "
                    f"collections"
                )
                exit(0)


@cli.command()
@click.pass_context
def create_cxg_artifacts(ctx):
    """
    Create cxg artifacts for all datasets in the database based on their explorer_url
    DO NOT run/use once dataset updates have shipped -- the s3 location will no longer be
    based on the explorer_url in all cases.
    You must first SSH into the target deployment using `make db/tunnel` before running.
    You must first set DEPLOYMENT_STAGE as an env var before running
    To run
    ./scripts/cxg_admin.py --deployment prod create-cxg-artifacts
    """
    with db_session_manager() as session:
        click.confirm(
            f"Are you sure you want to run this script? It will delete all of the current cxg artifacts and create new "
            f"ones based on the explorer_url?",
            abort=True,
        )
        session.query(DbDatasetArtifact).filter(DbDatasetArtifact.filetype == DatasetArtifactFileType.CXG).delete()
        session.commit()
        datasets = session.query(DbDataset.id, DbDataset.explorer_url).all()
        for dataset in datasets:
            if dataset.explorer_url:
                object_key = dataset.explorer_url.split("/")[-2]
                s3_uri = f"s3://{cxg_bucket.name}/{object_key}/"
                print(dataset.explorer_url, s3_uri)
                DatasetAsset.create(
                    session,
                    dataset_id=dataset.id,
                    filename="explorer_cxg",
                    filetype=DatasetArtifactFileType.CXG,
                    type_enum=DatasetArtifactType.REMIX,
                    user_submitted=True,
                    s3_uri=s3_uri,
                )


def get_database_uri() -> str:
    uri = urlparse(CorporaDbConfig().database_uri)
    uri = uri._replace(netloc="@".join([uri[1].split("@")[0], "localhost:5432"]))
    return uri.geturl()


if __name__ == "__main__":
    cli(obj={})
