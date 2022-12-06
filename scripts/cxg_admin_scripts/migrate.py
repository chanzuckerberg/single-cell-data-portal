import uuid
import click
import datetime
import logging
import requests
import os
import sys

from sqlalchemy import schema




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

def migrate_redesign_read(ctx):

    collections = []
    collection_versions = []
    datasets = []
    dataset_versions = []
    artifacts = []

    def strip_prefixes(v):
        if v is None:
            return None
        v = str(v)
        if "." in v and v.split(".")[1].isupper():
            return v.split(".")[1]
        else :
            return v

    def strip_prefixes_dict(d):
        r = {}
        for k,v in d.items():
            r[k] = strip_prefixes(v)
        return r

    with db_session_manager() as session:
        for record in session.query(DbCollection):
            print(record.id, record.name, record.visibility, record.links)

            if record.visibility == CollectionVisibility.PRIVATE and record.revision_of is None:
                # New, unpublished collection
                version_id = str(uuid.uuid4())
                collection_id = record.id

                collection = {
                    "id": record.id,
                    "version_id": None, # not mapped yet, this is private
                    "originally_published_at": None, # Not yet published
                    "tombstoned": False,
                }

                collections.append(collection)

            elif record.visibility == CollectionVisibility.PRIVATE and record.revision_of is not None:
                # Private revision of an existing collection
                version_id = record.id
                collection_id = record.revision_of
            else:
                # Published collection
                version_id = str(uuid.uuid4())
                collection_id = record.id

                collection = {
                    "id": record.id,
                    "version_id": version_id,
                    "originally_published_at": record.published_at,
                    "tombstoned": False,
                }

                collections.append(collection)

            metadata = {
                "name": record.name,
                "description": record.description,
                "contact_name": record.contact_name,
                "contact_email": record.contact_email,
                "links": [{"type": strip_prefixes(link.link_type), "url": link.link_url, "name": link.link_name} for link in record.links]
            }

            dataset_ids = []
            for record_dataset in record.datasets:
                dataset_version_id = str(uuid.uuid4())
                dataset = {
                    "dataset_id": record_dataset.id,
                    "dataset_version_id": dataset_version_id,
                    "published_at": record_dataset.published_at,
                }

                dataset_metadata = {
                    "name": record_dataset.name,
                    "schema_version": record_dataset.schema_version,
                    "organism": record_dataset.organism,
                    "tissue": record_dataset.tissue,
                    "assay": record_dataset.assay,
                    "disease": record_dataset.disease,
                    "sex": record_dataset.sex,
                    "self_reported_ethnicity": record_dataset.self_reported_ethnicity,
                    "development_stage": record_dataset.development_stage,
                    "cell_type": record_dataset.cell_type,
                    "cell_count": record_dataset.cell_count,
                    "schema_version": record_dataset.schema_version,
                    "mean_genes_per_cell": record_dataset.mean_genes_per_cell,
                    "batch_condition": record_dataset.batch_condition,
                    "suspension_type": record_dataset.suspension_type,
                    "donor_id": record_dataset.donor_id,
                    "is_primary_data": record_dataset.is_primary_data,
                    "x_approximate_distribution": record_dataset.x_approximate_distribution,
                }

                if record_dataset.processing_status is not None:
                    status = record_dataset.processing_status.__dict__
                else:
                    status = {}

                artifact_ids = []
                for record_artifact in record_dataset.artifacts:
                    artifact = {
                        "id": record_artifact.id,
                        "type": strip_prefixes(record_artifact.filetype),
                        "uri": record_artifact.s3_uri,
                    }
                    artifact_ids.append(record_artifact.id)
                    artifacts.append(artifact)


                dataset_version = {
                    "version_id": dataset_version_id,
                    "dataset_id": record_dataset.id,
                    "collection_id": version_id,
                    "metadata": dataset_metadata,
                    "artifacts": artifact_ids,
                    "status": strip_prefixes_dict(status),
                }

                dataset_ids.append(dataset_version_id)
                datasets.append(dataset)
                dataset_versions.append(dataset_version)
            
            version = {
                "version_id": version_id,
                "collection_id": collection_id,
                "metadata": metadata,
                "owner": record.owner,
                "publisher_metadata": record.publisher_metadata,
                "published_at": record.published_at,
                "datasets": dataset_ids,
            }

            collection_versions.append(version)

    import json
    with open("migration/collections.json", "w") as f:
        json.dump(collections, f, default=str)

    with open("migration/collection_versions.json", "w") as f:
        json.dump(collection_versions, f, default=str)

    with open("migration/datasets.json", "w") as f:
        json.dump(datasets, f, default=str)

    with open("migration/dataset_versions.json", "w") as f:
        json.dump(dataset_versions, f, default=str)

    with open("migration/dataset_artifacts.json", "w") as f:
        json.dump(artifacts, f, default=str)


def migrate_redesign_write(ctx):
    import json
    with open("migration/collections.json", "r") as f:
        collections = json.load(f)

    with open("migration/collection_versions.json", "r") as f:
        collection_versions = json.load(f)

    with open("migration/datasets.json", "r") as f:
        datasets = json.load(f)

    with open("migration/dataset_versions.json", "r") as f:
        dataset_versions = json.load(f)

    with open("migration/dataset_artifacts.json", "r") as f:
        artifacts = json.load(f)

    from backend.layers.persistence.orm import CollectionVersion
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker, session as sql_session
    from sqlalchemy.orm import Session


    from backend.layers.persistence.orm import (
        Collection as CollectionRow,
        CollectionVersion as CollectionVersionRow,
        Dataset as DatasetRow,
        DatasetVersion as DatasetVersionRow,
        DatasetArtifact as DatasetArtifactRow
    )

    # from backend.layers.persistence.orm import metadata_obj

    # database_uri = "postgresql://postgres:secret@localhost"
    database_pass = os.getenv("PGPASSWORD")
    database_uri = f"postgresql://dataportal:{database_pass}@localhost//ebezzi-redesign-full"

    # Uncomment for local
    # database_uri = f"postgresql://postgres:secret@localhost"
    engine = create_engine(database_uri, connect_args={"connect_timeout": 5})

    # engine.execute(schema.CreateSchema('persistence_schema'))
    # metadata_obj.create_all(bind=engine)

    with Session(engine) as session:

        for collection in collections:
            canonical_collection_row = CollectionRow(
                id=collection["id"],
                version_id=collection["version_id"],
                originally_published_at=collection.get("originally_published_at"),
                tombstoned=collection["tombstoned"],
            )

            session.add(canonical_collection_row)
        session.commit()

        for version in collection_versions:
            collection_version_row = CollectionVersionRow(
                                                        collection_id=version["collection_id"],
                                                        version_id=version["version_id"],
                                                        owner=version["owner"],
                                                        metadata=version["metadata"],
                                                        publisher_metadata=version["publisher_metadata"],
                                                        published_at=version["published_at"],
                                                        datasets=version["datasets"],
                                                        created_at=version.get("created_at"),
                                                    )
            session.add(collection_version_row)
        session.commit()

        for dataset in datasets:
            dataset_row = DatasetRow(
                dataset_id=dataset["dataset_id"],
                dataset_version_id=dataset["dataset_version_id"],
                published_at=dataset.get("published_at"),
            )

            session.add(dataset_row)
        session.commit()

        for dataset_version in dataset_versions:
            dataset_version_row = DatasetVersionRow(
                version_id=dataset_version["version_id"],
                dataset_id=dataset_version["dataset_id"],
                collection_id=dataset_version["collection_id"],
                metadata=dataset_version["metadata"],
                artifacts=dataset_version["artifacts"],
                status=dataset_version["status"],
                created_at=dataset_version.get("created_at"),
            )

            session.add(dataset_version_row)
        session.commit()



        for artifact in artifacts:
            artifact_row = DatasetArtifactRow(
                id=artifact["id"],
                type=artifact["type"],
                uri=artifact["uri"],
            )

            session.add(artifact_row)
        session.commit()
