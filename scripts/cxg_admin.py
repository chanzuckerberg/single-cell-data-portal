#!/usr/bin/env python
import json
import logging
import os
import sys
from datetime import datetime, time

import click
import requests
from click import Context
from furl import furl
from requests import HTTPError

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.corpora_config import CorporaDbConfig
from backend.corpora.common.utils.json import CustomJSONEncoder
from backend.corpora.common.utils.db_session import db_session_manager, DBSessionMaker
from backend.corpora.common.corpora_orm import (
    CollectionLinkType,
    CollectionVisibility,
    DbCollection,
    DbCollectionLink,
    DbDataset,
    DatasetArtifactFileType,
    DbDatasetArtifact,
    ProcessingStatus,
)
from backend.corpora.common.entities import DatasetAsset
from backend.corpora.common.entities.dataset import Dataset
from backend.corpora.common.entities.collection import Collection
from backend.corpora.common.utils.s3_buckets import buckets

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
        logger.info("Cannot run this script for prod. Aborting.")
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

        logger.info("Deletions complete!")


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
    Retrieve user ids from auth0 before running or ping an engineer on the team to check the id of the owner in the
    database
    To run (from repo root)
    ./scripts/cxg_admin.py --deployment prod transfer-collections $CURR_OWNER_ID $NEW_OWNER_ID
    """

    with db_session_manager() as session:
        collections = session.query(DbCollection).filter(DbCollection.owner == curr_owner).all()
        new_owner_collections_count = len(session.query(DbCollection).filter(DbCollection.owner == new_owner).all())
        if len(collections):
            click.confirm(
                f"Are you sure you want to update the owner of {len(collections)} collection"
                f"{'s' if len(collections) > 1 else ''} from {curr_owner} to "
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
                    f"Updated owner of collection for {updated} collections. {new_owner} is now the owner of "
                    f"{[[x.name, x.id] for x in collections]}"
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
            "Are you sure you want to run this script? It will delete all of the current cxg artifacts and create new "
            "ones based on the explorer_url?",
            abort=True,
        )
        session.query(DbDatasetArtifact).filter(DbDatasetArtifact.filetype == DatasetArtifactFileType.CXG).delete()
        session.commit()
        datasets = session.query(DbDataset.id, DbDataset.explorer_url).all()
        for dataset in datasets:
            if dataset.explorer_url:
                object_key = dataset.explorer_url.split("/")[-2]
                s3_uri = f"s3://{buckets.explorer_bucket.name}/{object_key}/"
                click.echo(dataset.explorer_url, s3_uri)
                DatasetAsset.create(
                    session,
                    dataset_id=dataset.id,
                    filename="explorer_cxg",
                    filetype=DatasetArtifactFileType.CXG,
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
            "Are you sure you want to run this script? It will assign schema_version to all the datasets",
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
            "Are you sure you want to run this script? It will assign published_at to "
            "all of the existing collections and datasets",
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

        logger.info("----- Finished migrating published_at for collections! -----")

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

        logger.info("----- Finished migrating published_at for datasets! -----")


@cli.command()
@click.pass_context
def populate_revised_at(ctx):
    """
    Populates `revised_at` for each existing collection and dataset with the
    current datetime (UTC). This is a one-off procedure since revised_at will
    be set for collections and datasets when they are updated.
    """

    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to run this script? It will assign revised_at to "
            "all of the existing collections and datasets",
            abort=True,
        )

        now = datetime.utcnow()

        # Collections
        for record in session.query(DbCollection):
            collection_id = record.id

            # Skip private collection, since revised_at will be populated on
            # publish if there is a change to the collection.
            if record.visibility == CollectionVisibility.PRIVATE:
                logger.info(f"SKIPPING - Collection is PRIVATE | collection.id: {collection_id}")
                continue

            logger.info(f"Setting revised_at for collection {collection_id}")
            record.revised_at = now

        logger.info("----- Finished populating revised_at for collections! -----")

        # Datasets
        for record in session.query(DbDataset):
            dataset_id = record.id

            # Skip private dataset, since revised_at will be populated on
            # publish if there are any changes.
            if record.collection_visibility == CollectionVisibility.PRIVATE:
                logger.info(f"SKIPPING - Dataset's parent collection is PRIVATE | dataset.id: {dataset_id}")
                continue

            logger.info(f"Setting revised_at for dataset {dataset_id}")
            record.revised_at = now

        logger.info("----- Finished populating revised_at for datasets! -----")


@cli.command()
@click.pass_context
def strip_all_collection_fields(ctx):
    """
    Strip all the `collection` string fields, so whitespace at the beginning and the end are removed.
    """
    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to run this script? It will strip whitespaces to all the string"
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
def backfill_processing_status_for_datasets(ctx):
    """
    Backfills the `dataset_processing_status` table for datasets that do not have a matching record.
    """
    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to run this script? It will assign dataset_processing_status "
            "to all datasets that are missing it",
            abort=True,
        )

        for record in session.query(DbDataset):
            dataset_id = record.id
            if record.processing_status.processing_status is None:
                record.processing_status.processing_status = ProcessingStatus.SUCCESS
                logger.warning(f"Setting processing status for dataset {dataset_id} {record.collection_id}")
            else:
                logger.warning(f"{dataset_id} processing status is fine")


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
            "Are you sure you want to run this script? It will add a trailing slash to "
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


@cli.command()
@click.argument("dataset_uuid")
@click.pass_context
def reprocess_seurat(ctx: Context, dataset_uuid: str) -> None:
    """
    Reconverts the specified dataset to Seurat format in place.
    :param ctx: command context
    :param dataset_uuid: UUID of dataset to reconvert to Seurat format
    """
    import boto3
    from time import time

    deployment = ctx.obj["deployment"]

    click.confirm(
        f"Are you sure you want to run this script? "
        f"It will reconvert and replace the dataset {dataset_uuid} to Seurat in the {deployment} environment.",
        abort=True,
    )

    aws_account_id = get_aws_account_id()
    deployment = ctx.obj["deployment"]
    happy_stack_name = get_happy_stack_name(deployment)

    payload = {"dataset_uuid": dataset_uuid}

    client = boto3.client("stepfunctions")
    response = client.start_execution(
        stateMachineArn=f"arn:aws:states:us-west-2:{aws_account_id}:stateMachine:dp-{happy_stack_name}-seurat-sfn",
        name=f"{dataset_uuid}-{int(time())}",
        input=json.dumps(payload),
    )

    click.echo(
        f"Step function executing: "
        f"https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/executions/details/"
        f"{response['executionArn']}"
    )


def get_aws_account_id() -> str:
    import boto3

    sts = boto3.client("sts")
    return sts.get_caller_identity()["Account"]


def get_happy_stack_name(deployment) -> str:
    """
    Returns the name of the Happy stack for the specified deployment
    Note: This will only work with deployment={dev,stage,prod} and will not work with rdev!
    :param deployment: dev, stage or prod
    :return:
    """
    return f"{deployment}-{deployment}stack"


auth0_apis = {
    "staging": "https://czi-cellxgene-dev.us.auth0.com",
    "dev": "https://czi-cellxgene-dev.us.auth0.com",
    "prod": "https://corpora-prod.auth0.com",
}


@cli.command()
@click.argument("access_token")
@click.pass_context
def update_curator_names(ctx, access_token):
    """Add the curator name to all collection based on the owner of the collection.

    ACCESS_TOKEN: Retrieved from Auth0 console or generated using the Client ID and Client Secret.
    The application must be authorized to access to the Auth0 Management API with the following permissions read:users
    read:user_idp_tokens.
    """

    auth0_api = auth0_apis[ctx.obj["deployment"]]

    def get_collections_without_curator():
        logger.info("Gathering collections with no curator.")
        with db_session_manager() as session:
            from sqlalchemy import null

            _filter = (Collection.table.curator_name == null()) | (Collection.table.curator_name == "")
            _owners = [
                result.owner for result in session.query(Collection.table.owner).filter(_filter).distinct().all()
            ]
        return _owners

    def get_owner_info_from_auth0(owner, access_token):
        logger.info(f"Retrieving curator info for owner:{owner} from Auth0.")
        url_manager = furl(
            f"{auth0_api}/api/v2/users/{owner}",
            query_params=dict(fields=",".join(["given_name", "family_name", "name"]), include_fields="true"),
        )
        response = requests.get(url_manager.url, headers={"Authorization": f"Bearer {access_token}"})
        response.raise_for_status()
        if response.headers["X-RateLimit-Remaining"] == 0:
            wait = time.time() - response.headers["X-RateLimit-Reset"]
            time.sleep(wait)
        body = response.json()
        name = f"{body.get('given_name', '')} {body.get('family_name', '')}".strip() or body.get("name", "")
        return name

    def update_database_curator_name(owner_id, owner_name):
        logger.info(f"Updating collections owned by {owner_id} with curator_name:{owner_name}.")
        with db_session_manager() as session:
            collections = session.query(DbCollection).filter(DbCollection.owner == owner_id).all()
            for collection in collections:
                collection.curator_name = owner_name
                if collection.visibility == CollectionVisibility.PUBLIC:
                    collection.data_submission_policy_version = "2.0"

    unique_owners = get_collections_without_curator()
    owners = dict()
    bad_owners = []
    for owner in unique_owners:
        try:
            owner_name = get_owner_info_from_auth0(owner, access_token)
            if owner_name:
                owners[owner] = owner_name
        except HTTPError:
            bad_owners.append(owner)
            logger.exception(f"Failed to fetch Auth0 info for owner:{owner}")
    for owner_id, owner_name in owners.items():
        update_database_curator_name(owner_id, owner_name)


@cli.command()
@click.pass_context
def add_publisher_metadata(ctx):
    """Add publisher metadata to the current records"""

    from backend.corpora.common.providers import crossref_provider

    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to run this script? It will populate publisher_metadata for all "
            "datasets. This will also do N calls to Crossref.",
            abort=True,
        )

        import time
        import traceback

        provider = crossref_provider.CrossrefProvider()

        for record in session.query(DbCollection):
            collection = Collection.get_collection(session, record.id, record.visibility)
            if not collection:
                continue
            collection_id = record.id
            dois = [link.link_url for link in collection.links if link.link_type == CollectionLinkType.DOI]
            normalized_doi = collection.get_doi()

            if record.publisher_metadata:
                print("Already has metadata, skipping...")
                continue

            if normalized_doi:
                if "/" not in normalized_doi:
                    continue

                try:
                    metadata = provider.fetch_metadata(normalized_doi)
                    record.publisher_metadata = metadata
                    print(collection_id, dois, normalized_doi)
                except Exception as e:
                    print(record.id, normalized_doi, record.name, e)
                    print(traceback.format_exc())
                    record.publisher_metadata = {"error": str(e.__cause__)}

                time.sleep(2)


@cli.command()
@click.pass_context
def wmg_get_s3_uris(ctx):
    """
    Print dataset ids and s3_uris for all datasets meeting the wmg criteria
    ./scripts/cxg_admin.py --deployment dev wmg-get-s3-uris
    """

    from backend.wmg.data.extract import get_dataset_s3_uris
    with db_session_manager() as session:
        s3_uris = get_dataset_s3_uris()
        print(s3_uris)


@cli.command()
@click.pass_context
def refresh_preprint_doi(ctx):
    """Add publisher metadata to the current records"""

    from backend.corpora.common.providers import crossref_provider

    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to run this script? It will populate publisher_metadata for all "
            "datasets. This will also do N calls to Crossref.",
            abort=True,
        )

        provider = crossref_provider.CrossrefProvider()

        for record in session.query(DbCollection):
            collection = Collection.get_collection(session, record.id, record.visibility)
            if not collection:
                continue
            collection_id = record.id
            normalized_doi = collection.get_doi()

            if record.publisher_metadata and record.publisher_metadata.get("is_preprint"):
                print(f"{normalized_doi}: is_preprint")

                try:
                    published_doi = provider.fetch_preprint_published_doi(normalized_doi)
                    print(collection_id, normalized_doi, published_doi)
                    doi_record = session.query(DbCollectionLink).filter(
                        DbCollectionLink.collection_id == record.id,
                        DbCollectionLink.collection_visibility == record.visibility,
                        DbCollectionLink.link_type == CollectionLinkType.DOI,
                    )
                    doi_record.link_url = published_doi
                except Exception as e:
                    print(e)


@cli.command()
@click.pass_context
def cxg_remaster(ctx):
    """Cxg remaster"""

    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to remaster all the cxgs?",
            abort=True,
        )

        import boto3

        client = boto3.client("stepfunctions")
        from time import time, sleep

        for record in session.query(DbDataset):
            if not record.tombstone:
                dataset = Dataset.get(session, dataset_uuid=record.id)
                artifacts = [a.s3_uri for a in dataset.artifacts if a.filetype == DatasetArtifactFileType.CXG]
                if len(artifacts) > 0:
                    cxg = artifacts[0]

                    p = urlparse(cxg)
                    bucket = p.hostname
                    dataset_id = p.path.strip("/").strip(".cxg")

                    if dataset_id == "2e5273bd-aa36-4478-8f6f-62fa0abcea43":
                        continue

                    if bucket != "hosted-cellxgene-dev":
                        continue

                    print(bucket, dataset_id)

                    input = {"dataset_uuid": dataset_id}

                    aws_account_id = get_aws_account_id()
                    deployment = ctx.obj["deployment"]
                    happy_stack_name = get_happy_stack_name(deployment)

                    response = client.start_execution(
                        stateMachineArn=f"arn:aws:states:us-west-2:{aws_account_id}:stateMachine:dp-"
                                        f"{happy_stack_name}-cxg-remaster-sfn",
                        name=f"{dataset_id}-{int(time())}",
                        input=json.dumps(input),
                    )

                    print(response["executionArn"])
                    sleep(1)


def get_database_uri() -> str:
    uri = urlparse(CorporaDbConfig().database_uri)
    uri = uri._replace(netloc="@".join([uri[1].split("@")[0], "localhost:5432"]))
    return uri.geturl()


if __name__ == "__main__":
    cli(obj={})
