import click
from click import Context
import json
import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.common.entities.collection import Collection
from backend.common.entities.dataset import Dataset

from backend.common.utils.db_session import db_session_manager
from backend.common.utils.json import CustomJSONEncoder


def tombstone_collection(ctx: Context, uuid: str):
    """
    Tombstones the collection specified by ID.

    Before running, create a tunnel to the database, e.g.:

        AWS_PROFILE=single-cell-prod DEPLOYMENT_STAGE=prod make db/tunnel

    Then run as:

        ./scripts/cxg_admin.py --deployment prod tombstone-collection 7edef704-f63a-462c-8636-4bc86a9472bd

    :param ctx: command context
    :param uuid: ID that identifies the collection to tombstone
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
