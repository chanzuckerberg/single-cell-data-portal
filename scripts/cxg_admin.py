#!/usr/bin/env python
import json
import logging
import os
import sys

import click
from click import Context

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

import requests as rq

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
            click.echo(
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
@click.argument("collection_name")
@click.pass_context
def delete_collections(ctx, collection_name):
    """
    Delete collections from data portal staging or dev by collection name.

    You must first SSH into the target deployment using `make db/tunnel` before running.
    You must first set DEPLOYMENT_STAGE as an env var before running
    To run
    ./scripts/cxg_admin.py --deployment dev delete-collections <collection_name>

    Examples of valid collection_name:
        - String with no spaces: ThisCollection
        - String with spaces: "This Collection"
    """

    if ctx.obj["deployment"] == "prod":
        logger.info(f"Cannot run this script for prod. Aborting.")
        exit(0)

    click.confirm(
        f"Are you sure you want to run this script? It will delete all of the "
        f"collections with the name '{collection_name}' from the {ctx.obj['deployment']} environment.",
        abort=True,
    )

    with db_session_manager() as session:
        collections = session.query(DbCollection).filter_by(name=collection_name).all()

        if not collections:
            logger.info(f"There are no collections with the name '{collection_name}'. Aborting.")
            exit(0)

        logger.info(f"There are {len(collections)} collections with the name '{collection_name}'")

        for c in collections:
            collection = Collection.get_collection(session, c.id, CollectionVisibility.PUBLIC, include_tombstones=True)
            if not collection:
                collection = Collection.get_collection(
                    session, c.id, CollectionVisibility.PRIVATE, include_tombstones=True
                )

            # Delete collection
            logger.info(f"Starting deletion of collection | name: {collection_name} | id: {c.id}")
            collection.delete()

        logger.info(f"Deletions complete!")


@cli.command()
@click.argument("uuid")
@click.pass_context
def tombstone_collection(ctx: Context, uuid: str):
    """
    Tombstones the collection specified by UUID.

    Before running, create a tunnel to the database, e.g.:

        AWS_PROFILE=single-cell-prod DEPLOYMENT_STAGE=prod make db/tunnel

    Then run as:

        ./scripts/cxg_admin.py --deployment prod tombstone-collection 7edef704-f63a-462c-8636-4bc86a9472bd

    :param ctx: command context
    :param uuid: UUID that identifies the collection to tombstone
    """

    with db_session_manager() as session:
        collection = Collection.get_collection(session, uuid, include_tombstones=True)

        if not collection:
            click.echo(f"Collection:{uuid} not found!")
            exit(0)

        click.echo(
            json.dumps(
                collection.to_dict(remove_attr=["datasets", "links", "genesets"]),
                sort_keys=True,
                indent=2,
                cls=CustomJSONEncoder,
            )
        )
        click.confirm(
            f"Are you sure you want to tombstone the collection:{uuid} from cellxgene:{ctx.obj['deployment']}?",
            abort=True,
        )

        collection.tombstone_collection()

        tombstoned = Collection.get_collection(session, uuid, include_tombstones=True)
        if tombstoned.tombstone:
            click.echo(f"Tombstoned collection:{uuid}")
            exit(0)
        else:
            click.echo(f"Failed to tombstone collection:{uuid}")
            exit(1)


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
            click.echo(
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
                click.echo(dataset.explorer_url, s3_uri)
                DatasetAsset.create(
                    session,
                    dataset_id=dataset.id,
                    filename="explorer_cxg",
                    filetype=DatasetArtifactFileType.CXG,
                    type_enum=DatasetArtifactType.REMIX,
                    user_submitted=True,
                    s3_uri=s3_uri,
                )


@cli.command()
@click.pass_context
def migrate_schema_version(ctx):
    """
    Populates `schema_version` for each existing dataset. Since the schema version only exists
    in the cxg file and we don't want to open them, we will call the cellxgene explorer endpoint
    which contains the version. This is a one-off procedure since new datasets will have
    the version already set.
    """

    with db_session_manager() as session:
        click.confirm(
            f"Are you sure you want to run this script? It will assign schema_version to all the datasets",
            abort=True,
        )
        for record in session.query(DbDataset):
            dataset_id = record.id
            explorer_url = urlparse(record.explorer_url)
            url = f"https://api.{explorer_url.netloc}/cellxgene{explorer_url.path}api/v0.2/config"
            res = rq.get(url).json()
            version = res["config"]["corpora_props"]["version"]["corpora_schema_version"]
            logger.info(f"Setting version for dataset {dataset_id} to {version}")
            record.schema_version = version


@cli.command()
@click.pass_context
def migrate_published_at(ctx):
    """
    Populates `published_at` for each existing collection and dataset. This is a
    one-off procedure since published_at will be set for collections and new
    datasets when they are first published.
    """

    with db_session_manager() as session:
        click.confirm(
            f"Are you sure you want to run this script? It will assign published_at to "
            f"all of the existing collections and datasets",
            abort=True,
        )
        # Collections
        for record in session.query(DbCollection):
            collection_id = record.id

            # Skip private collection, since published_at will be populated when published.
            if record.visibility == CollectionVisibility.PRIVATE:
                logger.info(f"SKIPPING - Collection is PRIVATE | collection.id: {collection_id}")
                continue

            # Skip if published_at already populated.
            if record.published_at is not None:
                logger.info(f"SKIPPING - Collection already has published_at | collection.id: {collection_id}")
                continue

            logger.info(f"Setting published_at for collection {collection_id}")
            collection_created_at = record.created_at
            record.published_at = collection_created_at

        logger.info(f"----- Finished migrating published_at for collections! -----")

        # Datasets
        for record in session.query(DbDataset):
            dataset_id = record.id

            # Skip private dataset, since published_at will be populated when published.
            if record.collection_visibility == CollectionVisibility.PRIVATE:
                logger.info(f"SKIPPING - Dataset's parent collection is PRIVATE | dataset.id: {dataset_id}")
                continue

            # Skip if published_at already populated.
            if record.published_at is not None:
                logger.info(f"SKIPPING - Dataset already has published_at | dataset.id: {dataset_id}")
                continue

            logger.info(f"Setting published_at for dataset {dataset_id}")
            dataset_created_at = record.created_at
            record.published_at = dataset_created_at

        logger.info(f"----- Finished migrating published_at for datasets! -----")


@cli.command()
@click.pass_context
def strip_all_collection_fields(ctx):
    """
    Strip all the `collection` string fields, so whitespace at the beginning and the end are removed.
    """
    with db_session_manager() as session:
        click.confirm(
            f"Are you sure you want to run this script? It will strip whitespaces to all the string"
            "fields in the `collection` table",
            abort=True,
        )

        query = """
            update project 
            set 
            owner=TRIM(owner), 
            name=TRIM(name), 
            description=TRIM(description), 
            contact_name=TRIM(contact_name), 
            contact_email=TRIM(contact_email), 
            data_submission_policy_version=TRIM(data_submission_policy_version)
            where
            owner != TRIM(owner) OR
            name != TRIM(name) OR
            description != TRIM(description) OR
            contact_name != TRIM(contact_name) OR
            contact_email != TRIM(contact_email) OR
            data_submission_policy_version != TRIM(data_submission_policy_version)
            """

        session.execute(query)
        session.commit()


@cli.command()
@click.pass_context
def add_trailing_slash_to_explorer_urls(ctx):
    """
    The explorer_url for datasets must end with a trailing slash to function
    properly. This script adds a trailing slash to a dataset's explorer_url
    if it already does not end with one.
    """
    with db_session_manager() as session:
        click.confirm(
            f"Are you sure you want to run this script? It will add a trailing slash to "
            "a dataset's explorer_url if it already does not end with one.",
            abort=True,
        )

        for record in session.query(DbDataset):
            dataset_id = record.id

            if record.explorer_url is None:
                logger.info(f"SKIPPING - Dataset does not have an explorer_url | dataset_id {dataset_id}")
                continue

            explorer_url = record.explorer_url.strip()
            if explorer_url[-1] == "/":
                logger.info(
                    f"SKIPPING - Dataset explorer_url already ends with trailing slash | dataset_id {dataset_id}"
                )
                continue

            logger.info(
                f"Adding trailing slash to dataset explorer_url | dataset_id {dataset_id} | "
                f"original explorer_url: {explorer_url}"
            )
            explorer_url_w_slash = explorer_url + "/"
            record.explorer_url = explorer_url_w_slash

        logger.info("----- Finished adding trailing slash to explorer_url for datasets! ----")


def get_database_uri() -> str:
    uri = urlparse(CorporaDbConfig().database_uri)
    uri = uri._replace(netloc="@".join([uri[1].split("@")[0], "localhost:5432"]))
    return uri.geturl()


if __name__ == "__main__":
    cli(obj={})
