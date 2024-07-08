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
from backend.common.providers.crossref_provider import CrossrefProvider
from backend.common.utils.aws import AwsSecret
from backend.layers.business.business import BusinessLogic
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.batch_job_provider import BatchJobProvider
from backend.layers.thirdparty.cloudfront_provider import CloudfrontProvider
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.step_function_provider import StepFunctionProvider
from backend.layers.thirdparty.uri_provider import UriProvider
from scripts.cxg_admin_scripts import (
    dataset_details,
    deletions,
    reprocess_datafile,
    request_logs,
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
    os.environ["DEPLOYMENT_STAGE"] = deployment
    happy_env = "stage" if deployment == "staging" else deployment
    stackname = f"{happy_env}stack"
    happy_config = json.loads(AwsSecret(f"happy/env-{happy_env}-config").value)
    os.environ["DATASETS_BUCKET"] = happy_config["s3_buckets"]["datasets"]["name"]

    ctx.obj["deployment"] = deployment
    ctx.obj["business_logic"] = BusinessLogic(
        DatabaseProvider(get_database_uri()),
        BatchJobProvider(),
        CrossrefProvider(),
        StepFunctionProvider(),
        S3Provider(),
        UriProvider(),
    )
    ctx.obj["stackname"] = stackname
    ctx.obj["cloudfront_provider"] = CloudfrontProvider()


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
@click.argument("collection_id")
@click.pass_context
def tombstone_collection(ctx: click.Context, collection_id: str):
    """
    Tombstones a public Collection specified by collection_id.
    To run:
        ./scripts/cxg_admin.py --deployment prod tombstone-collection 01234567-89ab-cdef-0123-456789abcdef

    :param ctx: command context
    :param collection_id: uuid that identifies the Collection to tombstone
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)  # Suppress type-related warnings from db operations
        tombstones.tombstone_collection(ctx, collection_id)


@cli.command()
@click.argument("collection_id")
@click.pass_context
def resurrect_collection(ctx: click.Context, collection_id: str):
    """
    Resurrects a tombstoned Collection specified by collection_id.
    To run:
        ./scripts/cxg_admin.py --deployment prod resurrect-collection 01234567-89ab-cdef-0123-456789abcdef

    :param ctx: command context
    :param collection_id: uuid that identifies the Collection to resurrect
    """
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


@cli.group("schema-migration")
@click.pass_context
def schema_migration_cli(ctx):
    """
    Commands for schema migration
    """
    deployment = ctx.obj["deployment"]
    happy_env = "stage" if deployment == "staging" else ctx.obj["deployment"]
    happy_config = json.loads(AwsSecret(f"happy/env-{happy_env}-config").value)
    os.environ["ARTIFACT_BUCKET"] = happy_config["s3_buckets"]["artifact"]["name"]


@schema_migration_cli.command()
@click.pass_context
@click.argument("report_path", type=click.Path(exists=True))
def rollback_datasets(ctx, report_path: str):
    """
    Used to rollback a datasets to a previous version.

    ./scripts/cxg_admin.py schema-migration --deployment dev rollback-dataset report.json
    """
    schema_migration.rollback_dataset(ctx, report_path)

    
@schema_migration_cli.command()
@click.pass_context
@click.argument("execution_id")
@click.argument("output_path", type=click.Path(writable=True), default=".")
def generate_report(ctx, execution_id: str, output_path: str):
    """
    Generates a report for the schema migration process.
    ./scripts/cxg_admin.py --deployment dev schema-migration generate-report execution_id
    """
    schema_migration.generate_report(ctx, execution_id, output_path, os.environ["ARTIFACT_BUCKET"])

    
@cli.command()
@click.pass_context
@click.argument("request_id")
@click.option("--hours", default=6, help="Number of hours to look back in the logs")
def get_request_logs(ctx, request_id: str, hours: int):
    """
    Get the requests from from AWS cloudwatch given the request_id

    ./scripts/cxg_admin.py --deployment dev get-request-logs <request_id>
    """
    print(json.dumps(request_logs.get(request_id, hours, ctx.obj["deployment"], ctx.obj["stackname"]), indent=4))



if __name__ == "__main__":
    cli(obj={})
