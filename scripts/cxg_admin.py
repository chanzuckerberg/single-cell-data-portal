#!/usr/bin/env python
import json
import logging
import os
import sys
import warnings

import click

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from urllib.parse import urlparse

from backend.common.corpora_config import CorporaDbConfig
from backend.common.utils.aws import AwsSecret
from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.crossref_provider import CrossrefProvider
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.step_function_provider import StepFunctionProvider
from backend.layers.thirdparty.uri_provider import UriProvider
from scripts.cxg_admin_scripts import (
    dataset_details,
    deletions,
    migrate,
    reprocess_datafile,
    schema_migration,
    tombstones,
    updates,
)

logging.basicConfig()
logger = logging.getLogger(__name__)

os.environ["CORPORA_LOCAL_DEV"] = "1"


def get_database_uri() -> str:
    uri = urlparse(CorporaDbConfig().database_uri)
    uri = uri._replace(netloc="@".join([uri[1].split("@")[0], "localhost:5432"]))
    return uri.geturl()


@click.group()
@click.option(
    "--deployment",
    default=lambda: os.environ.get("DEPLOYMENT_STAGE", "test"),
    show_default=True,
    help="The name of the deployment to target.",
)
@click.pass_context
def cli(ctx, deployment):
    """
    For all cxg_admin scripts
    You must first SSH into the target deployment using `make db/tunnel` before running (see
    backend/database readme for specific instructions)
    You must first set DEPLOYMENT_STAGE as an env var before running

    """
    if deployment == "test":
        return
    if deployment not in ("dev", "staging", "prod"):
        logging.error("The deployment arg must be one of 'dev', 'staging', or 'prod'")
        exit(1)
    happy_env = "stage" if deployment == "staging" else deployment
    happy_config = json.loads(AwsSecret(f"happy/env-{happy_env}-config").value)
    os.environ["DATASETS_BUCKET"] = happy_config["s3_buckets"]["datasets"]["name"]

    os.environ["DEPLOYMENT_STAGE"] = deployment
    ctx.obj["deployment"] = deployment
    ctx.obj["business_logic"] = BusinessLogic(
        DatabaseProvider(get_database_uri()), CrossrefProvider(), StepFunctionProvider(), S3Provider(), UriProvider()
    )


# Commands to delete artifacts (collections or datasets)


@cli.command()
@click.argument("id")
@click.pass_context
def delete_dataset(ctx, id):
    """Delete a dataset from Cellxgene."""
    deletions.delete_dataset(ctx, id)


@cli.command()
@click.argument("collection_name")
@click.pass_context
def delete_collections(ctx, collection_name):
    """
    Delete collections from data portal staging or dev by collection name.
    To run:
    ./scripts/cxg_admin.py --deployment dev delete-collections <collection_name>

    Examples of valid collection_name:
        - String with no spaces: ThisCollection
        - String with spaces: "This Collection"
    """
    deletions.delete_collections(ctx, collection_name)


# Commands to tombstone artifacts (datasets or collections)


@cli.command()
@click.argument("id")
@click.pass_context
def tombstone_collection(ctx: click.Context, id: str):
    """
    Tombstones the collection specified by ID.
    To run:
        ./scripts/cxg_admin.py --deployment prod tombstone-collection 7edef704-f63a-462c-8636-4bc86a9472bd

    :param ctx: command context
    :param id: ID that identifies the collection to tombstone
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # Suppress type-related warnings from db operations
        tombstones.tombstone_collection(ctx, id)


@cli.command()
@click.argument("collection_id")
@click.pass_context
def resurrect_collection(ctx: click.Context, collection_id: str):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # Suppress type-related warnings from db operations
        tombstones.resurrect_collection(ctx, collection_id)


@cli.command()
@click.argument("id")
@click.pass_context
def tombstone_dataset(ctx, id):
    """
    Remove a dataset from Cellxgene. This will delete its artifacts/genesets and mark the dataset as tombstoned so
     it no longer shows up in the data portal.
    To run:
      ./scripts/cxg_admin.py --deployment staging tombstone-dataset "57cf1b53-af10-49e5-9a86-4bc70d0c92b6"
    """
    tombstones.tombstone_dataset(ctx, id)


# Command to update different metadata fields


@cli.command()
@click.argument("collection_id")
@click.argument("new_owner")
@click.pass_context
def update_collection_owner(ctx, collection_id, new_owner):
    """Update the owner of a cellxgene collection.
    To run:
    ./scripts/cxg_admin.py --deployment prod update-collection-owner "$COLLECTION_ID $NEW_OWNER_ID
    """
    updates.update_collection_owner(ctx, collection_id, new_owner)


@cli.command()
@click.argument("curr_owner")
@click.argument("new_owner")
@click.pass_context
def transfer_collections(ctx, curr_owner, new_owner):
    """Transfer all collections owned by the curr_owner to the new_owner.
    Retrieve user ids from auth0 before running or ping an engineer on the team to check the id of the owner in the
    database
    To run:
    ./scripts/cxg_admin.py --deployment prod transfer-collections $CURR_OWNER_ID $NEW_OWNER_ID
    """
    updates.transfer_collections(ctx, curr_owner, new_owner)


@cli.command()
@click.pass_context
def strip_all_collection_fields(ctx):
    """
    Strip all the `collection` string fields, so whitespace at the beginning and the end are removed.
    """
    updates.strip_all_collection_fields(ctx)


@cli.command()
@click.pass_context
def add_trailing_slash_to_explorer_urls(ctx):
    """
    The explorer_url for datasets must end with a trailing slash to function
    properly. This script adds a trailing slash to a dataset's explorer_url
    if it already does not end with one.
    """
    updates.add_trailing_slash_to_explorer_urls(ctx)


@cli.command()
@click.argument("access_token")
@click.pass_context
def update_curator_names(ctx, access_token):
    """Add the curator name to all collection based on the owner of the collection.

    ACCESS_TOKEN: Retrieved from Auth0 console or generated using the Client ID and Client Secret.
    The application must be authorized to access to the Auth0 Management API with the following permissions read:users
    read:user_idp_tokens.
    """

    updates.update_curator_names(ctx, access_token)


@cli.command()
@click.pass_context
def add_publisher_metadata(ctx):
    """Add publisher metadata to the current records"""
    updates.add_publisher_metadata(ctx)


@cli.command()
@click.pass_context
def refresh_preprint_doi(ctx):
    """Add publisher metadata to the current records"""
    updates.refresh_preprint_doi(ctx)


# Commands to migrate the data, typically one off scripts run to populate db for existing rows after adding a new field


@cli.command()
@click.pass_context
def create_cxg_artifacts(ctx):
    """
    Create cxg artifacts for all datasets in the database based on their explorer_url
    DO NOT run/use once dataset updates have shipped -- the s3 location will no longer be
    based on the explorer_url in all cases.
    To run
    ./scripts/cxg_admin.py --deployment prod create-cxg-artifacts
    """
    migrate.create_cxg_artifacts(ctx)


@cli.command()
@click.pass_context
def migrate_schema_version(ctx):
    """
    Populates `schema_version` for each existing dataset. Since the schema version only exists
    in the cxg file and we don't want to open them, we will call the cellxgene explorer endpoint
    which contains the version. This is a one-off procedure since new datasets will have
    the version already set.
    """

    migrate.migrate_schema_version(ctx)


@cli.command()
@click.pass_context
def migrate_published_at(ctx):
    """
    Populates `published_at` for each existing collection and dataset. This is a
    one-off procedure since published_at will be set for collections and new
    datasets when they are first published.
    """
    migrate.migrate_published_at(ctx)


@cli.command()
@click.pass_context
def populate_revised_at(ctx):
    """
    Populates `revised_at` for each existing collection and dataset with the
    current datetime (UTC). This is a one-off procedure since revised_at will
    be set for collections and datasets when they are updated.
    """
    migrate.populate_revised_at(ctx)


@cli.command()
@click.pass_context
def backfill_processing_status_for_datasets(ctx):
    """
    Backfills the `dataset_processing_status` table for datasets that do not have a matching record.
    """
    migrate.backfill_processing_status_for_datasets(ctx)


# Commands to reprocess dataset artifacts (seurat or cxg)


@cli.command()
@click.argument("dataset_id")
@click.pass_context
def reprocess_seurat(ctx: click.Context, dataset_id: str) -> None:
    """
    Reconverts the specified dataset to Seurat format in place.
    :param ctx: command context
    :param dataset_id: ID of dataset to reconvert to Seurat format
    """
    reprocess_datafile.reprocess_seurat(ctx, dataset_id)


@cli.command()
@click.pass_context
def cxg_remaster(ctx):
    """Cxg remaster v2"""
    reprocess_datafile.cxg_remaster(ctx)


# Command to pull information from the db
@cli.command()
@click.pass_context
def wmg_get_asset_urls(ctx):
    """
    Print dataset ids and s3_uris for all datasets meeting the wmg criteria
    ./scripts/cxg_admin.py --deployment dev wmg-get-s3-uris
    """

    from backend.wmg.pipeline.integrated_corpus.extract import get_dataset_asset_urls

    asset_urls = get_dataset_asset_urls()
    print(asset_urls)


# Command to pull information from the db
@cli.command()
@click.pass_context
def get_public_datasets(ctx):
    """
    Print id, name, organism, tissue, assay, sex, cell_count, explorer_url, and S3 uris for all public datasets
    ./scripts/cxg_admin.py --deployment dev get-public-datasets
    """
    published_datasets = dataset_details.get_public_dataset_details()
    print(json.dumps(published_datasets, indent=2))


@cli.command()
@click.pass_context
def migrate_redesign_read(ctx):
    """
    Dumps the existing database to a set of .json files that match to the new schema. Files will be in the migration
    folder
    ./scripts/cxg_admin.py --deployment dev migrate-redesign-read
    """
    migrate.migrate_redesign_read(ctx)


@cli.command()
@click.pass_context
def migrate_redesign_write(ctx):
    """
    Reads the files generated by `migrate-redesign-read` and saves them to the desired environment
    ./scripts/cxg_admin.py --deployment dev migrate-redesign-write
    """
    migrate.migrate_redesign_write(ctx)


@cli.command()
@click.pass_context
def migrate_redesign_debug(ctx):
    """
    Used to debug the schema migration process. To be removed.
    ./scripts/cxg_admin.py --deployment dev migrate-redesign-debug
    """
    migrate.migrate_redesign_debug(ctx)


@cli.command()
@click.pass_context
def migrate_redesign_correct_published_at(ctx):
    """
    Used to debug the schema migration process. To be removed.
    ./scripts/cxg_admin.py --deployment dev migrate-redesign-debug
    """
    migrate.migrate_redesign_correct_published_at(ctx)


@cli.command()
@click.pass_context
@click.argument("report_patj", type=click.Path(exists=True))
def rollback_datasets(ctx, report_path: str):
    """
    Used to rollback a datasets to a previous version.

    ./scripts/cxg_admin.py --deployment dev rollback-dataset report.json
    """
    schema_migration.rollback_dataset(ctx, report_path)


if __name__ == "__main__":
    cli(obj={})
