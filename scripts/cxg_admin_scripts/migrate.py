import click
import datetime
import logging
import requests
import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.common.utils.db_session import db_session_manager
from urllib.parse import urlparse
from backend.common.corpora_orm import (
    CollectionVisibility,
    DbCollection,
    DbDataset,
    DatasetArtifactFileType,
    DbDatasetArtifact,
    ProcessingStatus,
)
from backend.common.entities import DatasetAsset
from backend.common import buckets

logging.basicConfig()
logger = logging.getLogger(__name__)


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
            res = requests.get(url).json()
            version = res["config"]["corpora_props"]["version"]["corpora_schema_version"]
            logger.info(f"Setting version for dataset {dataset_id} to {version}")
            record.schema_version = version


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
            if record.collection.visibility == CollectionVisibility.PRIVATE:
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
            if record.collection.visibility == CollectionVisibility.PRIVATE:
                logger.info(f"SKIPPING - Dataset's parent collection is PRIVATE | dataset.id: {dataset_id}")
                continue

            logger.info(f"Setting revised_at for dataset {dataset_id}")
            record.revised_at = now

        logger.info("----- Finished populating revised_at for datasets! -----")


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
