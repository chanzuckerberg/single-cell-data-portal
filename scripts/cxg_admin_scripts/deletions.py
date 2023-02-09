import json
import logging
import os
import sys

import click

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.common.corpora_orm import CollectionVisibility, DbCollection
from backend.common.entities.collection import Collection
from backend.common.entities.dataset import Dataset
from backend.common.utils.db_session import db_session_manager
from backend.common.utils.json import CustomJSONEncoder

logging.basicConfig()
logger = logging.getLogger(__name__)


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
