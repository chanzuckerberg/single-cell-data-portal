import uuid
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

def migrate_redesign(ctx):

    collections = []
    collection_versions = []
    datasets = []
    dataset_versions = []
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
                "links": [{"type": link.link_type, "url": link.link_url, "name": link.link_name} for link in record.links]
            }

            dataset_ids = []
            for dataset in record.datasets:
                pass

            # class Dataset(CanonicalDataset):

#     __table__ = Table(
#         "Dataset",
#         mapper_registry.metadata,
#         Column("dataset_id", Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)),
#         Column("dataset_version_id", Column(UUID(as_uuid=True), default=uuid.uuid4)),
#         Column("published_at", Column(DateTime))
#     )


# @mapper_registry.mapped
# class DatasetVersion(DatasetVersionModel):

#     artifacts: List[DatasetArtifactId] = field(default=list())
#     canonical_dataset: CanonicalDataset = field(default=None)

#     __table__ = Table(
#         "DatasetVersion",
#         mapper_registry.metadata,
#         Column("version_id", Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)),
#         Column("dataset_id", Column(UUID(as_uuid=True), default=uuid.uuid4)),
#         Column("collection_id", Column(UUID(as_uuid=True), default=uuid.uuid4)),
#         Column("metadata", Column(JSON)),
#         Column("artifacts", Column(ARRAY(UUID(as_uuid=True)))),
#         Column("status", Column(JSON))
#     )


# @mapper_registry.mapped
# class DatasetArtifact(DatasetArtifactModel):

#     __table__ = Table(
#         "DatasetArtifact",
#         mapper_registry.metadata,
#         Column("id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
#         Column("type", Enum(DatasetArtifactType)),
#         Column("uri", String)
#     )


            version = {
                "version_id": version_id,
                "collection_id": collection_id,
                "metadata": metadata,
                "owner": record.owner,
                "publisher_metadata": record.publisher_metadata,
                "published_at": record.published_at,
                "datasets": [],
            }


            collection_versions.append(version)

    import json
    with open("migration/collections.json", "w") as f:
        json.dump(collections, f, default=str)

    with open("migration/collection_versions.json", "w") as f:
        json.dump(collection_versions, f, default=str)



            # for dataset in record.datasets:
                # print(dataset.id)


# id                             | 43d4bb39-21af-4d05-b973-4c1fed7b916c
# owner                          | google-oauth2|106687116196602342983
# name                           | Transcriptional Programming of Normal and Inflamed Human Epidermis at Single-Cell Resolution
# description                    | Perturbations in the transcriptional programs specifying epidermal differentiation cause diverse skin pathologies ranging from impaired barrier function to inflammatory skin disease. However, the global scope and organization of this complex cellular program remain undefined. Here we report single-cell RNA sequencing profiles of 92,889 human epidermal cells from 9 normal and 3 inflamed skin samples. Transcriptomics-derived keratinocyte subpopulations reflect classic epidermal strata but also sharply compartmentalize epithelial functions such as cell-cell communication, inflammation, and WNT pathway modulation. In keratinocytes, ∼12% of assessed transcript expression varies in coordinate patterns, revealing undescribed gene expression programs governing epidermal homeostasis. We also identify molecular fingerprints of inflammatory skin states, including S100 activation in the interfollicular epidermis of normal scalp, enrichment of a CD1C+CD301A+ myeloid dendritic cell population in psoriatic epidermis, and IL1βhiCCL3hiCD14+ monocyte-derived macrophages enriched in foreskin. This compendium of RNA profiles provides a critical step toward elucidating epidermal diseases of development, differentiation, and inflammation.
# created_at                     | 2022-10-04 19:52:10.600886
# updated_at                     | 2022-10-04 19:52:10.600891
# visibility                     | PRIVATE
# contact_email                  | raymond.cho@ucsf.edu
# contact_name                   | Raymond J. Cho
# data_submission_policy_version |
# tombstone                      | f
# published_at                   |
# revised_at                     |
# curator_name                   | Batuhan Cakir
# publisher_metadata             | {"authors": [{"given": "Jeffrey B.", "family": "Cheng"}, {"given": "Andrew J.", "family": "Sedgewick"}, {"given": "Alex I.", "family": "Finnegan"}, {"given": "Paymann", "family": "Harirchian"}, {"given": "Jerry", "family": "Lee"}, {"given": "Sunjong", "family": "Kwon"}, {"given": "Marlys S.", "family": "Fassett"}, {"given": "Justin", "family": "Golovato"}, {"given": "Matthew", "family": "Gray"}, {"given": "Ruby", "family": "Ghadially"}, {"given": "Wilson", "family": "Liao"}, {"given": "Bethany E.", "family": "Perez White"}, {"given": "Theodora M.", "family": "Mauro"}, {"given": "Thaddeus", "family": "Mully"}, {"given": "Esther A.", "family": "Kim"}, {"given": "Hani", "family": "Sbitany"}, {"given": "Isaac M.", "family": "Neuhaus"}, {"given": "Roy C.", "family": "Grekin"}, {"given": "Siegrid S.", "family": "Yu"}, {"given": "Joe W.", "family": "Gray"}, {"given": "Elizabeth", "family": "Purdom"}, {"given": "Ralf", "family": "Paus"}, {"given": "Charles J.", "family": "Vaske"}, {"given": "Stephen C.", "family": "Benz"}, {"given": "Jun S.", "family": "Song"}, {"given": "Raymond J.", "family": "Cho"}], "published_year": 2018, "published_month": 10, "published_day": 1, "published_at": 1538352000.0, "journal": "Cell Reports", "is_preprint": false}
# revision_of                    |



        # for record in session.query()
            
            
            # CollectionMetadata(
            #     name=record.name,
            #     description=record.description
            # )





#     __table__ = Table(
#         "Collection",
#         mapper_registry.metadata,
#         Column("id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
#         Column("version_id", UUID(as_uuid=True), default=uuid.uuid4),
#         Column("originally_published_at", Column(DateTime)),
#         Column("tombstoned", Column(BOOLEAN))
#     )


# @mapper_registry.mapped
# class CollectionVersion(CollectionVersionModel):

#     canonical_collection: CanonicalCollection = field(default=None)

#     __table__ = Table(
#         "CollectionVersion",
#         mapper_registry.metadata,
#         Column("version_id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4),
#         Column("collection_id", UUID(as_uuid=True), default=uuid.uuid4),
#         Column("metadata", Column(JSON)),
#         Column("owner", Column(String)),
#         Column("publisher_metadata", Column(JSON)),
#         Column("published_at", Column(DateTime)),
#         Column("datasets", Column(ARRAY(UUID(as_uuid=True))))
#     )