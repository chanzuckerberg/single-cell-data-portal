import os
import sys

from click import Context

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa


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
    print(f"Testing DJH context: {ctx} uuid: {uuid}")
    import pprint

    pprint.pprint(ctx.obj)
    print(os.environ["DEPLOYMENT_STAGE"])
    pass


def tombstone_dataset(ctx, uuid):
    """
    Remove a dataset from Cellxgene. This will delete its artifacts/genesets and mark the dataset as tombstoned so
     it no longer shows up in the data portal.
     You must first SSH into the target deployment using `make db/tunnel` before running.
      ./scripts/cxg_admin.py --deployment staging tombstone-dataset "57cf1b53-af10-49e5-9a86-4bc70d0c92b6"

    """
    pass
